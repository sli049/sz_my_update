/**
 *  @file sz.c
 *  @author Sheng Di and Dingwen Tao
 *  @date Aug, 2016
 *  @brief SZ_Init, Compression and Decompression functions
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "sz.h"
#include "CompressElement.h"
#include "DynamicByteArray.h"
#include "DynamicIntArray.h"
#include "TightDataPointStorageD.h"
#include "TightDataPointStorageF.h"
#include "zlib.h"
#include "rw.h"
#include "Huffman.h"
#include "conf.h"
#include "utility.h"
//#include "CurveFillingCompressStorage.h"

int versionNumber[4] = {SZ_VER_MAJOR,SZ_VER_MINOR,SZ_VER_BUILD,SZ_VER_REVISION};
//int SZ_SIZE_TYPE = 8;


float delta_t_opt[5] = {3.4630957287598745e-05, 5.2311700705000656e-05, 3.5391612072822954e-05, 5.3441990932428505e-05, 5.4088709466407214e-05};
//sihuan added, this is only to test ideas. need to be calculated online when ideas are implemented

int dataEndianType = LITTLE_ENDIAN_DATA; //*endian type of the data read from disk
int sysEndianType; //*sysEndianType is actually set automatically.

//the confparams should be separate between compression and decopmression, in case of mutual-affection when calling compression/decompression alternatively
sz_params *confparams_cpr = NULL; //used for compression
sz_params *confparams_dec = NULL; //used for decompression 

sz_exedata *exe_params = NULL;

int sz_with_regression = SZ_WITH_LINEAR_REGRESSION; //SZ_NO_REGRESSION

/*following global variables are desgined for time-series based compression*/
/*sz_varset is not used in the single-snapshot data compression*/
SZ_VarSet* sz_varset = NULL;
sz_multisteps *multisteps = NULL;
sz_tsc_metadata *sz_tsc = NULL;
SZ_Variable* v_global = NULL; //sihuan make it global
int TheCurVarNum = 0;//sihuan make a global var counter

//only for Pastri compressor
#ifdef PASTRI
pastri_params pastri_par;
#endif

HuffmanTree* SZ_Reset()
{
	return createDefaultHuffmanTree();
}

int SZ_Init(const char *configFilePath)
{
	int loadFileResult = SZ_LoadConf(configFilePath);
	if(loadFileResult==SZ_NSCS)
		return SZ_NSCS;
	
	exe_params->SZ_SIZE_TYPE = sizeof(size_t);
	
	if(confparams_cpr->szMode == SZ_TEMPORAL_COMPRESSION)
	{
		initSZ_TSC();
	}
	return SZ_SCES;
}

int SZ_Init_Params(sz_params *params)
{
	int x = 1;
	char *y = (char*)&x;
	int endianType = BIG_ENDIAN_SYSTEM;
	if(*y==1) endianType = LITTLE_ENDIAN_SYSTEM;

	sysEndianType = endianType;
	exe_params->SZ_SIZE_TYPE = sizeof(size_t);

	// set default values
	if(params->max_quant_intervals > 0) 
		params->maxRangeRadius = params->max_quant_intervals/2;
	else
		params->max_quant_intervals = params->maxRangeRadius*2;

	exe_params->intvCapacity = params->maxRangeRadius*2;
	exe_params->intvRadius = params->maxRangeRadius;

	if(params->quantization_intervals>0)
	{
		updateQuantizationInfo(params->quantization_intervals);
		exe_params->optQuantMode = 0;
	}
	else
		exe_params->optQuantMode = 1;


	if(params->quantization_intervals%2!=0)
	{
		printf("Error: quantization_intervals must be an even number!\n");
		return SZ_NSCS;
	}

	confparams_cpr = (sz_params*)malloc(sizeof(sz_params));
	memcpy(confparams_cpr, params, sizeof(sz_params));	

	return SZ_SCES;
}

int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	int dimension;
	if(r1==0) 
	{
		dimension = 0;
	}
	else if(r2==0) 
	{
		dimension = 1;
	}
	else if(r3==0) 
	{
		dimension = 2;
	}
	else if(r4==0) 
	{
		dimension = 3;
	}
	else if(r5==0) 
	{
		dimension = 4;
	}
	else 
	{
		dimension = 5;
	}
	return dimension;	
}

size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	size_t dataLength;
	if(r1==0) 
	{
		dataLength = 0;
	}
	else if(r2==0) 
	{
		dataLength = r1;
	}
	else if(r3==0) 
	{
		dataLength = r1*r2;
	}
	else if(r4==0) 
	{
		dataLength = r1*r2*r3;
	}
	else if(r5==0) 
	{
		dataLength = r1*r2*r3*r4;
	}
	else 
	{
		dataLength = r1*r2*r3*r4*r5;
	}
	return dataLength;
}

/*-------------------------------------------------------------------------*/
/**
    @brief      Perform Compression 
    @param      data           data to be compressed
    @param      outSize        the size (in bytes) after compression
    @param		r5,r4,r3,r2,r1	the sizes of each dimension (supporting only 5 dimensions at most in this version.
    @return     compressed data (in binary stream) or NULL(0) if any errors

 **/
/*-------------------------------------------------------------------------*/
unsigned char* SZ_compress_args(int dataType, void *data, size_t *outSize, int errBoundMode, double absErrBound, 
double relBoundRatio, double pwrBoundRatio, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	//TODO
	confparams_cpr->dataType = dataType;
	if(dataType==SZ_FLOAT)
	{
		unsigned char *newByteData = NULL;
		
		SZ_compress_args_float(&newByteData, (float *)data, r5, r4, r3, r2, r1, 
		outSize, errBoundMode, absErrBound, relBoundRatio, pwrBoundRatio);
		
		return newByteData;
	}
	else if(dataType==SZ_DOUBLE)
	{
		unsigned char *newByteData;
		SZ_compress_args_double(&newByteData, (double *)data, r5, r4, r3, r2, r1, 
		outSize, errBoundMode, absErrBound, relBoundRatio, pwrBoundRatio);
		
		return newByteData;
	}
	else if(dataType==SZ_INT64)
	{
		unsigned char *newByteData;
		SZ_compress_args_int64(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;
	}		
	else if(dataType==SZ_INT32) //int type
	{
		unsigned char *newByteData;
		SZ_compress_args_int32(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;
	}
	else if(dataType==SZ_INT16)
	{
		unsigned char *newByteData;
		SZ_compress_args_int16(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;		
	}
	else if(dataType==SZ_INT8)
	{
		unsigned char *newByteData;
		SZ_compress_args_int8(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;
	}
	else if(dataType==SZ_UINT64)
	{
		unsigned char *newByteData;
		SZ_compress_args_uint64(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;
	}		
	else if(dataType==SZ_UINT32) //int type
	{
		unsigned char *newByteData;
		SZ_compress_args_uint32(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;
	}
	else if(dataType==SZ_UINT16)
	{
		unsigned char *newByteData;
		SZ_compress_args_uint16(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;		
	}
	else if(dataType==SZ_UINT8)
	{
		unsigned char *newByteData;
		SZ_compress_args_uint8(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;
	} 	
	else
	{
		printf("Error: dataType can only be SZ_FLOAT, SZ_DOUBLE, SZ_INT8/16/32/64 or SZ_UINT8/16/32/64.\n");
		return NULL;
	}
}

int SZ_compress_args2(int dataType, void *data, unsigned char* compressed_bytes, size_t *outSize, 
int errBoundMode, double absErrBound, double relBoundRatio, double pwrBoundRatio, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	unsigned char* bytes = SZ_compress_args(dataType, data, outSize, errBoundMode, absErrBound, relBoundRatio, pwrBoundRatio, r5, r4, r3, r2, r1);
    memcpy(compressed_bytes, bytes, *outSize);
    free(bytes); 
	return SZ_SCES;
}

int SZ_compress_args3(int dataType, void *data, unsigned char* compressed_bytes, size_t *outSize, int errBoundMode, double absErrBound, double relBoundRatio, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1,
size_t s5, size_t s4, size_t s3, size_t s2, size_t s1,
size_t e5, size_t e4, size_t e3, size_t e2, size_t e1)
{
	confparams_cpr->dataType = dataType;
	if(dataType==SZ_FLOAT)
	{
		SZ_compress_args_float_subblock(compressed_bytes, (float *)data, 
		r5, r4, r3, r2, r1,
		s5, s4, s3, s2, s1,
		e5, e4, e3, e2, e1,
		outSize, errBoundMode, absErrBound, relBoundRatio);
		
		return SZ_SCES;
	}
	else if(dataType==SZ_DOUBLE)
	{
		SZ_compress_args_double_subblock(compressed_bytes, (double *)data, 
		r5, r4, r3, r2, r1,
		s5, s4, s3, s2, s1,
		e5, e4, e3, e2, e1,
		outSize, errBoundMode, absErrBound, relBoundRatio);
		
		return SZ_SCES;
	}
	else
	{
		printf("Error (in SZ_compress_args3): dataType can only be SZ_FLOAT or SZ_DOUBLE.\n");
		return SZ_NSCS;
	}	
}

unsigned char *SZ_compress(int dataType, void *data, size_t *outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{	
	unsigned char *newByteData = SZ_compress_args(dataType, data, outSize, confparams_cpr->errorBoundMode, confparams_cpr->absErrBound, confparams_cpr->relBoundRatio, 
	confparams_cpr->pw_relBoundRatio, r5, r4, r3, r2, r1);
	return newByteData;
}

//////////////////
/*-------------------------------------------------------------------------*/
/**
    @brief      Perform Compression 
    @param      data           data to be compressed
    @param		reservedValue  the reserved value
    @param      outSize        the size (in bytes) after compression
    @param		r5,r4,r3,r2,r1	the sizes of each dimension (supporting only 5 dimensions at most in this version.
    @return     compressed data (in binary stream)

 **/
/*-------------------------------------------------------------------------*/
unsigned char *SZ_compress_rev_args(int dataType, void *data, void *reservedValue, size_t *outSize, int errBoundMode, double absErrBound, double relBoundRatio, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	unsigned char *newByteData;
	//TODO
	printf("SZ compression with reserved data is TO BE DONE LATER.\n");
	exit(0);
	
	return newByteData;	
}

int SZ_compress_rev_args2(int dataType, void *data, void *reservedValue, unsigned char* compressed_bytes, size_t *outSize, int errBoundMode, double absErrBound, double relBoundRatio, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	confparams_cpr->dataType = dataType;
	unsigned char* bytes = SZ_compress_rev_args(dataType, data, reservedValue, outSize, errBoundMode, absErrBound, relBoundRatio, r5, r4, r3, r2, r1);
	memcpy(compressed_bytes, bytes, *outSize);
	free(bytes); //free(bytes) is removed , because of dump error at MIRA system (PPC architecture), fixed?
	return 0;
}

unsigned char *SZ_compress_rev(int dataType, void *data, void *reservedValue, size_t *outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	unsigned char *newByteData;
	//TODO
	printf("SZ compression with reserved data is TO BE DONE LATER.\n");
	exit(0);
	
	return newByteData;
}

void *SZ_decompress(int dataType, unsigned char *bytes, size_t byteLength, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	if(confparams_dec==NULL)
		confparams_dec = (sz_params*)malloc(sizeof(sz_params));
	memset(confparams_dec, 0, sizeof(sz_params));
	if(exe_params==NULL)
		exe_params = (sz_exedata*)malloc(sizeof(sz_exedata));
	memset(exe_params, 0, sizeof(sz_exedata));
	
	int x = 1;
	char *y = (char*)&x;
	if(*y==1)
		sysEndianType = LITTLE_ENDIAN_SYSTEM;
	else //=0
		sysEndianType = BIG_ENDIAN_SYSTEM;
	
	if(dataType == SZ_FLOAT)
	{
		float *newFloatData;
		SZ_decompress_args_float(&newFloatData, r5, r4, r3, r2, r1, bytes, byteLength);
		return newFloatData;	
	}
	else if(dataType == SZ_DOUBLE)
	{
		double *newDoubleData;
		SZ_decompress_args_double(&newDoubleData, r5, r4, r3, r2, r1, bytes, byteLength);
		return newDoubleData;	
	}
	else if(dataType == SZ_INT8)
	{
		int8_t *newInt8Data;
		SZ_decompress_args_int8(&newInt8Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newInt8Data;
	}
	else if(dataType == SZ_INT16)
	{
		int16_t *newInt16Data;
		SZ_decompress_args_int16(&newInt16Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newInt16Data;
	}
	else if(dataType == SZ_INT32)
	{
		int32_t *newInt32Data;
		SZ_decompress_args_int32(&newInt32Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newInt32Data;
	}
	else if(dataType == SZ_INT64)
	{
		int64_t *newInt64Data;
		SZ_decompress_args_int64(&newInt64Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newInt64Data;
	}
	else if(dataType == SZ_UINT8)
	{
		uint8_t *newUInt8Data;
		SZ_decompress_args_uint8(&newUInt8Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newUInt8Data;
	}
	else if(dataType == SZ_UINT16)
	{
		uint16_t *newUInt16Data;
		SZ_decompress_args_uint16(&newUInt16Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newUInt16Data;
	}
	else if(dataType == SZ_UINT32)
	{
		uint32_t *newUInt32Data;
		SZ_decompress_args_uint32(&newUInt32Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newUInt32Data;
	}
	else if(dataType == SZ_UINT64)
	{
		uint64_t *newUInt64Data;
		SZ_decompress_args_uint64(&newUInt64Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newUInt64Data;
	}
	else 
	{
		printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
		return NULL;	
	}
}

/**
 * 
 * 
 * return number of elements or -1 if any errors
 * */
size_t SZ_decompress_args(int dataType, unsigned char *bytes, size_t byteLength, void* decompressed_array, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	//size_t i;
	size_t nbEle = computeDataLength(r5,r4,r3,r2,r1);
	
	if(dataType == SZ_FLOAT)
	{
		float* data = (float *)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		float* data_array = (float *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(float));
		//for(i=0;i<nbEle;i++)
		//	data_array[i] = data[i];	
		free(data); //this free operation seems to not work with BlueG/Q system.	
	}
	else if (dataType == SZ_DOUBLE)
	{
		double* data = (double *)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		double* data_array = (double *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(double));
		//for(i=0;i<nbEle;i++)
		//	data_array[i] = data[i];
		free(data); //this free operation seems to not work with BlueG/Q system.	
	}
	else if(dataType == SZ_INT8)
	{
		int8_t* data = (int8_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		int8_t* data_array = (int8_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(int8_t));
		free(data);
	}
	else if(dataType == SZ_INT16)
	{
		int16_t* data = (int16_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		int16_t* data_array = (int16_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(int16_t));
		free(data);	
	}
	else if(dataType == SZ_INT32)
	{
		int32_t* data = (int32_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		int32_t* data_array = (int32_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(int32_t));
		free(data);	
	}
	else if(dataType == SZ_INT64)
	{
		int64_t* data = (int64_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		int64_t* data_array = (int64_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(int64_t));
		free(data);		
	}
	else if(dataType == SZ_UINT8)
	{
		uint8_t* data = (uint8_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		uint8_t* data_array = (uint8_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(uint8_t));
		free(data);
	}
	else if(dataType == SZ_UINT16)
	{
		uint16_t* data = (uint16_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		uint16_t* data_array = (uint16_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(uint16_t));
		free(data);		
	}
	else if(dataType == SZ_UINT32)
	{
		uint32_t* data = (uint32_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		uint32_t* data_array = (uint32_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(uint32_t));
		free(data);		
	}
	else if(dataType == SZ_UINT64)
	{
		uint64_t* data = (uint64_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		uint64_t* data_array = (uint64_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(uint64_t));
		free(data);			
	}
	else
	{ 
		printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
		return SZ_NSCS; //indicating error		
	}

	return nbEle;
}


sz_metadata* SZ_getMetadata(unsigned char* bytes)
{
	int index = 0, i, isConstant, isLossless;
	size_t dataSeriesLength = 0;
	int versions[3] = {0,0,0};
	for (i = 0; i < 3; i++)
		versions[i] = bytes[index++]; //3
	unsigned char sameRByte = bytes[index++]; //1
	isConstant = sameRByte & 0x01;
	//confparams_dec->szMode = (sameRByte & 0x06)>>1;
	isLossless = (sameRByte & 0x10)>>4;
	if(exe_params==NULL)
	{
		exe_params = (sz_exedata *)malloc(sizeof(struct sz_exedata));
		memset(exe_params, 0, sizeof(struct sz_exedata));
	}
	exe_params->SZ_SIZE_TYPE = ((sameRByte & 0x40)>>6)==1?8:4;
	
	sz_params* params = convertBytesToSZParams(&(bytes[index]));
	if(confparams_dec!=NULL)
		free(confparams_dec);
	confparams_dec = params;	
	index += MetaDataByteLength;
	
	if(params->dataType!=SZ_FLOAT && params->dataType!= SZ_DOUBLE) //if this type is an Int type
		index++; //jump to the dataLength info byte address
	dataSeriesLength = bytesToSize(&(bytes[index]));// 4 or 8	
	index += exe_params->SZ_SIZE_TYPE;
	index += 4; //max_quant_intervals
	
	sz_metadata* metadata = (sz_metadata*)malloc(sizeof(struct sz_metadata));
	
	metadata->versionNumber[0] = versions[0];
	metadata->versionNumber[1] = versions[1];
	metadata->versionNumber[2] = versions[2];
	metadata->isConstant = isConstant;
	metadata->isLossless = isLossless;
	metadata->sizeType = exe_params->SZ_SIZE_TYPE;
	metadata->dataSeriesLength = dataSeriesLength;
	
	metadata->conf_params = confparams_dec;
	
	int defactoNBBins = 0; //real # bins
	if(isConstant==0 && isLossless==0)
	{
		int radExpoL = 0, segmentL = 0, pwrErrBoundBytesL = 0;
		if(metadata->conf_params->errorBoundMode >= PW_REL)
		{
			radExpoL = 1;
			segmentL = exe_params->SZ_SIZE_TYPE;
			pwrErrBoundBytesL = 4;
		}
		
		int offset_typearray = 3 + 1 + MetaDataByteLength + exe_params->SZ_SIZE_TYPE + 4 + radExpoL + segmentL + pwrErrBoundBytesL + 4 + (4 + params->dataType*4) + 1 + 8 
				+ exe_params->SZ_SIZE_TYPE + exe_params->SZ_SIZE_TYPE + exe_params->SZ_SIZE_TYPE + 4;
		defactoNBBins = bytesToInt_bigEndian(bytes+offset_typearray);
	}
	
	metadata->defactoNBBins = defactoNBBins;
	return metadata;
}

void SZ_printMetadata(sz_metadata* metadata)
{
	printf("=================SZ Compression Meta Data=================\n");
	printf("Version:                        \t %d.%d.%d\n", metadata->versionNumber[0], metadata->versionNumber[1], metadata->versionNumber[2]);
	printf("Constant data?:                 \t %s\n", metadata->isConstant==1?"YES":"NO");
	printf("Lossless?:                      \t %s\n", metadata->isLossless==1?"YES":"NO");
	printf("Size type (size of # elements): \t %d bytes\n", metadata->sizeType); 
	printf("Num of elements:                \t %zu\n", metadata->dataSeriesLength);
		
	sz_params* params = metadata->conf_params;
	
	switch(params->dataType)
	{
	case SZ_FLOAT:
		printf("Data type:                      \t FLOAT\n");
		break;
	case SZ_DOUBLE:
		printf("Data type:                      \t DOUBLE\n");
		break;
	case SZ_INT8:
		printf("Data type:                      \t INT8\n");
		break;	
	case SZ_INT16:
		printf("Data type:                      \t INT16\n");
		break;
	case SZ_INT32:
		printf("Data type:                      \t INT32\n");
		break;	
	case SZ_INT64:
		printf("Data type:                      \t INT64\n");
		break;	
	case SZ_UINT8:
		printf("Data type:                      \t UINT8\n");
		break;	
	case SZ_UINT16:
		printf("Data type:                      \t UINT16\n");
		break;
	case SZ_UINT32:
		printf("Data type:                      \t UINT32\n");
		break;	
	case SZ_UINT64:
		printf("Data type:                      \t UINT64\n");
		break;				
	}
	
	if(exe_params->optQuantMode==1)
	{
		printf("quantization_intervals:         \t 0\n");
		printf("max_quant_intervals:            \t %d\n", params->max_quant_intervals);
		printf("actual used # intervals:        \t %d\n", metadata->defactoNBBins);
	}
	else
	{
		printf("quantization_intervals:         \t %d\n", params->quantization_intervals);
		printf("max_quant_intervals:            \t - %d\n", params->max_quant_intervals);		
	}
	
	printf("dataEndianType (prior raw data):\t %s\n", dataEndianType==BIG_ENDIAN_DATA?"BIG_ENDIAN":"LITTLE_ENDIAN");
	printf("sysEndianType (at compression): \t %s\n", sysEndianType==1?"BIG_ENDIAN":"LITTLE_ENDIAN");
	printf("sampleDistance:                 \t %d\n", params->sampleDistance);
	printf("predThreshold:                  \t %f\n", params->predThreshold);
	switch(params->szMode)
	{
	case SZ_BEST_SPEED:
		printf("szMode:                         \t SZ_BEST_SPEED (without Gzip)\n");
		break;
	case SZ_BEST_COMPRESSION:
		printf("szMode:                         \t SZ_BEST_COMPRESSION (with Gzip)\n");
		break;
	}
	switch(params->gzipMode)
	{
	case Z_BEST_SPEED:
		printf("gzipMode:                       \t Z_BEST_SPEED\n");
		break;
	case Z_DEFAULT_COMPRESSION:
		printf("gzipMode:                       \t Z_BEST_SPEED\n");
		break;	
	case Z_BEST_COMPRESSION:
		printf("gzipMode:                       \t Z_BEST_COMPRESSION\n");
		break;
	}
	
	switch(params->errorBoundMode)
	{
	case ABS:
		printf("errBoundMode:                   \t ABS\n");
		printf("absErrBound:                    \t %f\n", params->absErrBound);
		break;
	case REL:
		printf("errBoundMode:                   \t REL (based on value_range extent)\n");
		printf("relBoundRatio:                  \t %f\n", params->relBoundRatio);
		break;
	case ABS_AND_REL:
		printf("errBoundMode:                   \t ABS_AND_REL\n");
		printf("absErrBound:                    \t %f\n", params->absErrBound);
		printf("relBoundRatio:                  \t %f\n", params->relBoundRatio);
		break;
	case ABS_OR_REL:
		printf("errBoundMode:                   \t ABS_OR_REL\n");
		printf("absErrBound:                    \t %f\n", params->absErrBound);
		printf("relBoundRatio:                  \t %f\n", params->relBoundRatio);
		break;
	case PSNR:
		printf("errBoundMode:                   \t PSNR\n");
		printf("psnr:                           \t %f\n", params->psnr);
		break;
	case PW_REL:
		printf("errBoundMode:                   \t PW_REL\n");
		break;
	case ABS_AND_PW_REL:
		printf("errBoundMode:                   \t ABS_AND_PW_REL\n");
		printf("absErrBound:                    \t %f\n", params->absErrBound);
		break;
	case ABS_OR_PW_REL:
		printf("errBoundMode:                   \t ABS_OR_PW_REL\n");
		printf("absErrBound:                    \t %f\n", params->absErrBound);
		break;
	case REL_AND_PW_REL:
		printf("errBoundMode:                   \t REL_AND_PW_REL\n");
		printf("range_relBoundRatio:            \t %f\n", params->relBoundRatio);
		break;
	case REL_OR_PW_REL:
		printf("errBoundMode:                   \t REL_OR_PW_REL\n");
		printf("range_relBoundRatio:            \t %f\n", params->relBoundRatio);
		break;
	}
	
	if(params->errorBoundMode>=PW_REL && params->errorBoundMode<=REL_OR_PW_REL)
	{
		printf("pw_relBoundRatio:               \t %f\n", params->pw_relBoundRatio);
		printf("segment_size:                   \t %d\n", params->segment_size);
		switch(params->pwr_type)
		{
		case SZ_PWR_MIN_TYPE:
			printf("pwrType:                    \t SZ_PWR_MIN_TYPE\n");
			break;
		case SZ_PWR_AVG_TYPE:
			printf("pwrType:                    \t SZ_PWR_AVG_TYPE\n");
			break;
		case SZ_PWR_MAX_TYPE:
			printf("pwrType:                    \t SZ_PWR_MAX_TYPE\n");
			break;
		}
	}
}

/*-----------------------------------batch data compression--------------------------------------*/

void filloutDimArray(size_t* dim, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	if(r2==0)
		dim[0] = r1;
	else if(r3==0)
	{
		dim[0] = r2;
		dim[1] = r1;
	}
	else if(r4==0)
	{
		dim[0] = r3;
		dim[1] = r2;
		dim[2] = r1;
	}
	else if(r5==0)
	{
		dim[0] = r4;
		dim[1] = r3;
		dim[2] = r2;
		dim[3] = r1;
	}
	else
	{
		dim[0] = r5;
		dim[1] = r4;
		dim[2] = r3;
		dim[3] = r2;
		dim[4] = r1;		
	}
}

size_t compute_total_batch_size()
{
	size_t eleNum = 0, totalSize = 0;
	SZ_Variable* p = sz_varset->header;
	while(p->next!=NULL)
	{
		eleNum = computeDataLength(p->next->r5, p->next->r4, p->next->r3, p->next->r2, p->next->r1);
		if(p->next->dataType==SZ_FLOAT)
			totalSize += (eleNum*4);
		else
			totalSize += (eleNum*8);
		p=p->next;
	}
	return totalSize;
}

int isZlibFormat(unsigned char magic1, unsigned char magic2)
{
	if(magic1==104&&magic2==5) //DC+BS
		return 1;
	if(magic1==104&&magic2==129) //DC+DC
		return 1;
	if(magic1==104&&magic2==222) //DC+BC
		return 1;
	if(magic1==120&&magic2==1) //BC+BS
		return 1;
	if(magic1==120&&magic2==156) //BC+DC
		return 1;
	if(magic1==120&&magic2==218) //BC+BS
		return 1;
	return 0;
}

void SZ_registerVar(char* varName, int dataType, void* data, 
			int errBoundMode, double absErrBound, double relBoundRatio, double pwRelBoundRatio, 
			size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	if(sz_tsc==NULL)
		initSZ_TSC();
		
	char str[256];
	SZ_batchAddVar(varName, dataType, data, 
			errBoundMode, absErrBound, relBoundRatio, pwRelBoundRatio, r5, r4, r3, r2, r1);
	sprintf(str, "%d: %s : %zuX%zuX%zuX%zu%zu : %d : %f : %f : %f\n", sz_varset->count - 1, varName, r5, r4, r3, r2, r1, errBoundMode, absErrBound, relBoundRatio, pwRelBoundRatio);
	fputs(str, sz_tsc->metadata_file);
}

int SZ_deregisterVar(char* varName)
{
	int state = SZ_batchDelVar(varName);
	return state;
}

#ifdef HAVE_TIMECMPR
int SZ_compress_ts(unsigned char** newByteData, size_t *outSize)
{
	confparams_cpr->szMode = SZ_TEMPORAL_COMPRESSION;
	confparams_cpr->predictionMode = SZ_PREVIOUS_VALUE_ESTIMATE;
	
	SZ_VarSet* vset = sz_varset;
	size_t *outSize_ = (size_t*)malloc(sizeof(size_t)*vset->count);
	memset(outSize_, 0, sizeof(size_t)*vset->count);
	unsigned char** compressBuffer = (unsigned char**)malloc(vset->count*sizeof(unsigned char*));//to store compressed bytes
	
	char *metadata_str = (char*)malloc(vset->count*256);
	memset(metadata_str, 0, vset->count*256);
	sprintf(metadata_str, "step %d", sz_tsc->currentStep);
	
	int i = 0, totalSize = 0;
	SZ_Variable* v = NULL;
	reorder_vars(vset);//sihuan added

	size_t preLen = sz_tsc->intersect_size;
	v = vset->header->next;
	size_t dataLen = computeDataLength(v->r5, v->r4, v->r3, v->r2, v->r1);
	//printf("data length is: %zu\n", dataLen);
	size_t cur_intersect_size;

	//#if 0
	if (sz_tsc->currentStep % confparams_cpr->snapshotCmprStep == 0){
		cur_intersect_size = dataLen;
		//sihuan added: should write the reordered input for later decompression validation
		write_reordered_tofile(vset, dataLen);

		int64_t* tmp_hist_index = (int64_t*) malloc(sizeof(int64_t)*cur_intersect_size);
		v = vset->header->next;
		while(strcmp(v->varName, "index")) v = v->next;
		sz_tsc->hist_index = tmp_hist_index;
		memcpy(sz_tsc->hist_index, (int64_t*)v->data, sizeof(int64_t)*cur_intersect_size);
	}
	else {
		printf("prelen, dataLen before calling intersect: %zu, %zu\n", preLen, dataLen);
		unsigned char* tmp_bitarray = (unsigned char*) malloc(sizeof(unsigned char)*preLen);
		sz_tsc->bit_array = tmp_bitarray;
		cur_intersect_size = intersectAndsort(sz_tsc->hist_index, preLen, vset, dataLen, sz_tsc->bit_array);
		//sihuan added: should write the reordered input for later decompression validation
		write_reordered_tofile(vset, dataLen);
		//now write the bitarray to an output file
		char bitarr_out_name[50];
		sprintf(bitarr_out_name, "bitarray_out_sz_%d.sz1", sz_tsc->currentStep);
		unsigned char* bit_array_zlib_cmprs;
		unsigned long bit_array_zlib_cmprs_length = zlib_compress(sz_tsc->bit_array, preLen, &bit_array_zlib_cmprs, 1);
		unsigned char* bit_array_zlib_cmprs_meta = (unsigned char*) malloc(bit_array_zlib_cmprs_length + sizeof(size_t) + sizeof(unsigned long));
		//need also write the original size and the compressed byte length to the output file
		sizeToBytes(bit_array_zlib_cmprs_meta, preLen);
		//bit_array_zlib_cmprs_meta += sizeof(size_t);
		longToBytes_bigEndian(bit_array_zlib_cmprs_meta + sizeof(size_t), bit_array_zlib_cmprs_length);
		//bit_array_zlib_cmprs_meta += sizeof(unsigned long);
		memcpy(bit_array_zlib_cmprs_meta + sizeof(size_t) + sizeof(unsigned long), bit_array_zlib_cmprs, bit_array_zlib_cmprs_length);

		int status_tmp;
		writeByteData(bit_array_zlib_cmprs_meta, bit_array_zlib_cmprs_length+sizeof(size_t)+sizeof(unsigned long), bitarr_out_name, &status_tmp);
		free(bit_array_zlib_cmprs);
		free(bit_array_zlib_cmprs_meta);

		if ((sz_tsc->currentStep+1) % confparams_cpr->snapshotCmprStep == 0) free(sz_tsc->hist_index);
		else {
			free(sz_tsc->hist_index);
			int64_t* tmp_hist_index = (int64_t*) malloc(sizeof(int64_t)*cur_intersect_size);
			v = vset->header->next;
			while(strcmp(v->varName, "index")) v = v->next;
			sz_tsc->hist_index = tmp_hist_index;
			memcpy(sz_tsc->hist_index, (int64_t*)v->data, sizeof(int64_t)*cur_intersect_size);
		}
		//printf("sz_tsc->intersection_size is %zu", sz_tsc->intersect_size);
	}
	sz_tsc->intersect_size = cur_intersect_size;

	//#endif

	int partial = 0;

	if (sz_tsc->currentStep % confparams_cpr->snapshotCmprStep != 0) partial = 1; //time-based compression is partial

	for(i=0;i<vset->count-1;i++)//sihuan modified, the last variable is id which is not needed to compress. so that is why do count-1
	{
		if (i == 0) v = vset->header->next;
		multisteps = v->multisteps; //assign the v's multisteps to the global variable 'multisteps', which will be used in the following compression.
		//sihuan debug
		#if 0
		printf("current variable name: %s\n", v->varName);
		int tmp = 0;
		for (tmp = 0; tmp < 5; tmp++){
			if (strcmp(v->varName, "index")) printf("%.10f  ", ((float*)v->data)[tmp]);
			else printf("%lld  ", ((int64_t*)v->data)[tmp]);
		}
		printf("end printing the variable: %s\n", v->varName);
		#endif

		if(v->dataType==SZ_FLOAT)
		{
			if (partial){
				size_t outSize_phase1, outSize_phase2;
				unsigned char* buffer1 = NULL;
				unsigned char* buffer2 = NULL;
				SZ_compress_args_float_ps(&buffer1, (float*)v->data, v->r5, v->r4, v->r3, v->r2, v->r1, &outSize_phase1, v->errBoundMode, v->absErrBound, v->relBoundRatio, v->pwRelBoundRatio, partial, 1);
				SZ_compress_args_float_ps(&buffer2, (float*)v->data, v->r5, v->r4, v->r3, v->r2, v->r1, &outSize_phase2, v->errBoundMode, v->absErrBound, v->relBoundRatio, v->pwRelBoundRatio, partial, 2);
				
				outSize_[i] = sizeof(size_t) + outSize_phase1 + outSize_phase2;

				unsigned char* p = compressBuffer[i] = (unsigned char*)malloc(outSize_[i]);
				sizeToBytes(p, outSize_phase1);
				p += sizeof(size_t);
				memcpy(p, buffer1, outSize_phase1);
				p += outSize_phase1;
				memcpy(p, buffer2, outSize_phase2);//these 2 step memcpy may be further optimized
				//how to handle compressBuffer[i]?
				free(buffer1); free(buffer2);
				multisteps->compressionType = 2; //this type means the data is compressed both time-based and snapshot-based.
											  //the former part is done by time and the later is done by snapshot.
				//printf("The current varialbe of the current step %d is compressed by method: %d\n", sz_tsc->currentStep, multisteps->compressionType);

			}
			else{
				SZ_compress_args_float(&(compressBuffer[i]), (float*)v->data, v->r5, v->r4, v->r3, v->r2, v->r1, &outSize_[i], v->errBoundMode, v->absErrBound, v->relBoundRatio, v->pwRelBoundRatio);
				//printf("The current varialbe of the current step %d is compressed by method: %d\n", sz_tsc->currentStep, multisteps->compressionType);
			}
			//need to handle the rest part now

		}
		else if(v->dataType==SZ_DOUBLE)
		{
			SZ_compress_args_double(&(compressBuffer[i]), (double*)v->data, v->r5, v->r4, v->r3, v->r2, v->r1, &outSize_[i], v->errBoundMode, v->absErrBound, v->relBoundRatio, v->pwRelBoundRatio);
		}
		sprintf(metadata_str, "%s:%d,%d,%zu", metadata_str, i, multisteps->lastSnapshotStep, outSize_[i]);
		
		totalSize += outSize_[i];
		v->compressType = multisteps->compressionType;
		v = v->next;
	}
	
	sprintf(metadata_str, "%s\n", metadata_str);
	fputs(metadata_str, sz_tsc->metadata_file);
	free(metadata_str);
	
	//sizeof(int)==current time step; 2*sizeof(char)+sizeof(size_t)=={compressionType + datatype + compression_data_size}; 
	//sizeof(char)==# variables
	//sihuan added: sizeof(size_t)
	*outSize = sizeof(int) + sizeof(unsigned short) + sizeof(size_t) + totalSize + vset->count*(2*sizeof(unsigned char)+sizeof(size_t));
	*newByteData = (unsigned char*)malloc(*outSize); 
	unsigned char* p = *newByteData;

	intToBytes_bigEndian(p, sz_tsc->currentStep);
	p+=4;
	//shortToBytes(p, vset->count);
	shortToBytes(p, vset->count-1);// sihuan updated: subtract the id varialble: so vset->count-1;
	p+=2;
	sizeToBytes(p, sz_tsc->intersect_size);//sihuan added
	p += sizeof(size_t);//sihuan added
	
	//for(i=0;i<vset->count;i++)
	for(i=0;i<vset->count-1;i++)//sihuan updated: subtract the id variable: so vset->count-1;
	{
		//SZ_Variable* v = vset->header->next; //the code is correct? ed?
		if (i == 0) v = vset->header->next;
	
		*p = (unsigned char)v->compressType; //1 byte
		p++;
		*p = (unsigned char)v->dataType; //1 byte
		p++;
		sizeToBytes(p, outSize_[i]); //size_t
		p += sizeof(size_t);
		//sizeToBytes(p, v->r5); //size_t
		//p += sizeof(size_t);
		//sizeToBytes(p, v->r4); //size_t
		//p += sizeof(size_t);
		//sizeToBytes(p, v->r3); //size_t
		//p += sizeof(size_t);
		//sizeToBytes(p, v->r2); //size_t
		//p += sizeof(size_t);
		//sizeToBytes(p, v->r1); //size_t
		//p += sizeof(size_t);								
		memcpy(p, compressBuffer[i], outSize_[i]); //outSize_[i]
		p += outSize_[i];
		v = v->next; //sihuan corrected?
		//free(compressBuffer[i]); //do we need to free it here? why cannot free it?
	}

	//sihuan output size print
	printf("***********************printing results********************\n");
	for (i = 0; i < vset->count - 1; i++)
		printf("variable %d has compressed size: %zu, compression ratio: %.5f\n", i, outSize_[i], (float)(dataLen*sizeof(float))/(float)outSize_[i]);
	printf("***********************ending printing*********************\n");


	sz_tsc->currentStep ++;	
	free(outSize_);
	
	return SZ_SCES;
}

//#ifdef HAVE_TIMECMPR
int SZ_compress_ts_vlct(unsigned char** newByteData, size_t *outSize)
{
	TheCurVarNum = -1; //sihuan: finishing a snapshot should reset the var number.
	confparams_cpr->szMode = SZ_TEMPORAL_COMPRESSION;
	confparams_cpr->predictionMode = SZ_PREVIOUS_VALUE_ESTIMATE;
	
	SZ_VarSet* vset = sz_varset;
	size_t *outSize_ = (size_t*)malloc(sizeof(size_t)*vset->count);
	memset(outSize_, 0, sizeof(size_t)*vset->count);
	unsigned char** compressBuffer = (unsigned char**)malloc(vset->count*sizeof(unsigned char*));//to store compressed bytes
	
	char *metadata_str = (char*)malloc(vset->count*256);
	memset(metadata_str, 0, vset->count*256);
	sprintf(metadata_str, "step %d", sz_tsc->currentStep);
	
	int i = 0, totalSize = 0;
	SZ_Variable* v = NULL; //sihuan changed
	clock_t start, end;
	start = clock();
	reorder_vars(vset);//sihuan added
	end = clock();
	printf("sorting the id, swapping variables time: %.5f\n", (double)(end-start)/CLOCKS_PER_SEC);

	size_t preLen = sz_tsc->intersect_size;
	v = vset->header->next;
	size_t dataLen = computeDataLength(v->r5, v->r4, v->r3, v->r2, v->r1);
	//printf("data length is: %zu\n", dataLen);
	size_t cur_intersect_size;

	//#if 0
	printf("confparams_cpr->snapshotCmprStep is %d\n", confparams_cpr->snapshotCmprStep);
	if (sz_tsc->currentStep % confparams_cpr->snapshotCmprStep == 0){
		cur_intersect_size = dataLen;
		//sihuan added: should write the reordered input for later decompression validation
		write_reordered_tofile(vset, dataLen);

		int64_t* tmp_hist_index = (int64_t*) malloc(sizeof(int64_t)*cur_intersect_size);
		v = vset->header->next;
		while(strcmp(v->varName, "index")) v = v->next;
		sz_tsc->hist_index = tmp_hist_index;
		memcpy(sz_tsc->hist_index, (int64_t*)v->data, sizeof(int64_t)*cur_intersect_size);
	}
	else {
		printf("prelen, dataLen before calling intersect: %zu, %zu\n", preLen, dataLen);
		unsigned char* tmp_bitarray = (unsigned char*) malloc(sizeof(unsigned char)*preLen);
		sz_tsc->bit_array = tmp_bitarray;
		cur_intersect_size = intersectAndsort(sz_tsc->hist_index, preLen, vset, dataLen, sz_tsc->bit_array);
		//sihuan added: should calculate delta_t now;
		//delta_t_opt[sz_tsc->currentStep-1] = calculate_delta_t(cur_intersect_size);
		//sihuan added: should write the reordered input for later decompression validation
		write_reordered_tofile(vset, dataLen);
		//now write the bitarray to an output file
		char bitarr_out_name[50];
		sprintf(bitarr_out_name, "bitarray_out_sz_%d.sz1", sz_tsc->currentStep);
		unsigned char* bit_array_zlib_cmprs;
		unsigned long bit_array_zlib_cmprs_length = zlib_compress(sz_tsc->bit_array, preLen, &bit_array_zlib_cmprs, 1);
		unsigned char* bit_array_zlib_cmprs_meta = (unsigned char*) malloc(bit_array_zlib_cmprs_length + sizeof(size_t) + sizeof(unsigned long));
		//need also write the original size and the compressed byte length to the output file
		sizeToBytes(bit_array_zlib_cmprs_meta, preLen);
		//bit_array_zlib_cmprs_meta += sizeof(size_t);
		longToBytes_bigEndian(bit_array_zlib_cmprs_meta + sizeof(size_t), bit_array_zlib_cmprs_length);
		//bit_array_zlib_cmprs_meta += sizeof(unsigned long);
		memcpy(bit_array_zlib_cmprs_meta + sizeof(size_t) + sizeof(unsigned long), bit_array_zlib_cmprs, bit_array_zlib_cmprs_length);

		int status_tmp;
		writeByteData(bit_array_zlib_cmprs_meta, bit_array_zlib_cmprs_length+sizeof(size_t)+sizeof(unsigned long), bitarr_out_name, &status_tmp);
		free(bit_array_zlib_cmprs);
		free(bit_array_zlib_cmprs_meta);

		if ((sz_tsc->currentStep+1) % confparams_cpr->snapshotCmprStep == 0) free(sz_tsc->hist_index);
		else {
			free(sz_tsc->hist_index);
			int64_t* tmp_hist_index = (int64_t*) malloc(sizeof(int64_t)*cur_intersect_size);
			v = vset->header->next;
			while(strcmp(v->varName, "index")) v = v->next;
			sz_tsc->hist_index = tmp_hist_index;
			memcpy(sz_tsc->hist_index, (int64_t*)v->data, sizeof(int64_t)*cur_intersect_size);
		}
		//printf("sz_tsc->intersection_size is %zu", sz_tsc->intersect_size);
	}
	sz_tsc->intersect_size = cur_intersect_size;

	//#endif

	int partial = 0;

	if (sz_tsc->currentStep % confparams_cpr->snapshotCmprStep != 0) partial = 1; //time-based compression is partial

	//sihuan changed: need first find the variable "vx";
	//v = vset->header->next;
	//while(strcmp(v->varName, "vx")) v = v->next;


	for(i=0;i<vset->count-1;i++)//sihuan modified, the last variable is id which is not needed to compress. so that is why do count-1
	{
		TheCurVarNum++;//sihuan added
		if (i == 0) v = vset->header->next;
		v_global = v;
		multisteps = v->multisteps; //assign the v's multisteps to the global variable 'multisteps', which will be used in the following compression.
		//sihuan debug
		#if 0
		printf("current variable name: %s\n", v->varName);
		int tmp = 0;
		for (tmp = 0; tmp < 5; tmp++){
			if (strcmp(v->varName, "index")) printf("%.10f  ", ((float*)v->data)[tmp]);
			else printf("%lld  ", ((int64_t*)v->data)[tmp]);
		}
		printf("end printing the variable: %s\n", v->varName);
		#endif

		if(v->dataType==SZ_FLOAT)
		{
			if (partial){
				size_t outSize_phase1, outSize_phase2;
				unsigned char* buffer1 = NULL;
				unsigned char* buffer2 = NULL;
				SZ_compress_args_float_ps(&buffer1, (float*)v->data, v->r5, v->r4, v->r3, v->r2, v->r1, &outSize_phase1, v->errBoundMode, v->absErrBound, v->relBoundRatio, v->pwRelBoundRatio, partial, 1);
				SZ_compress_args_float_ps(&buffer2, (float*)v->data, v->r5, v->r4, v->r3, v->r2, v->r1, &outSize_phase2, v->errBoundMode, v->absErrBound, v->relBoundRatio, v->pwRelBoundRatio, partial, 2);
				
				outSize_[i] = sizeof(size_t) + outSize_phase1 + outSize_phase2;

				unsigned char* p = compressBuffer[i] = (unsigned char*)malloc(outSize_[i]);
				sizeToBytes(p, outSize_phase1);
				p += sizeof(size_t);
				memcpy(p, buffer1, outSize_phase1);
				p += outSize_phase1;
				memcpy(p, buffer2, outSize_phase2);//these 2 step memcpy may be further optimized
				//how to handle compressBuffer[i]?
				free(buffer1); free(buffer2);
				multisteps->compressionType = 2; //this type means the data is compressed both time-based and snapshot-based.
											  //the former part is done by time and the later is done by snapshot.
				//printf("The current varialbe of the current step %d is compressed by method: %d\n", sz_tsc->currentStep, multisteps->compressionType);

			}
			else{
				SZ_compress_args_float(&(compressBuffer[i]), (float*)v->data, v->r5, v->r4, v->r3, v->r2, v->r1, &outSize_[i], v->errBoundMode, v->absErrBound, v->relBoundRatio, v->pwRelBoundRatio);
				//printf("The current varialbe of the current step %d is compressed by method: %d\n", sz_tsc->currentStep, multisteps->compressionType);
			}
			//need to handle the rest part now

		}
		else if(v->dataType==SZ_DOUBLE)
		{
			SZ_compress_args_double(&(compressBuffer[i]), (double*)v->data, v->r5, v->r4, v->r3, v->r2, v->r1, &outSize_[i], v->errBoundMode, v->absErrBound, v->relBoundRatio, v->pwRelBoundRatio);
		}
		sprintf(metadata_str, "%s:%d,%d,%zu", metadata_str, i, multisteps->lastSnapshotStep, outSize_[i]);
		
		totalSize += outSize_[i];
		v->compressType = multisteps->compressionType;
		v = v->next;
	}

	//v = vset->header->next;
#if 0
	for(i=0;i<3;i++)//sihuan modified, the last variable is id which is not needed to compress. so that is why do count-1
	{
		TheCurVarNum++;//sihuan added
		//if (i == 0) v = vset->header->next;
		v_global = v;
		multisteps = v->multisteps; //assign the v's multisteps to the global variable 'multisteps', which will be used in the following compression.
		//sihuan debug
		#if 0
		printf("current variable name: %s\n", v->varName);
		int tmp = 0;
		for (tmp = 0; tmp < 5; tmp++){
			if (strcmp(v->varName, "index")) printf("%.10f  ", ((float*)v->data)[tmp]);
			else printf("%lld  ", ((int64_t*)v->data)[tmp]);
		}
		printf("end printing the variable: %s\n", v->varName);
		#endif

		if(v->dataType==SZ_FLOAT)
		{
			if (partial){
				size_t outSize_phase1, outSize_phase2;
				unsigned char* buffer1 = NULL;
				unsigned char* buffer2 = NULL;
				SZ_compress_args_float_ps(&buffer1, (float*)v->data, v->r5, v->r4, v->r3, v->r2, v->r1, &outSize_phase1, v->errBoundMode, v->absErrBound, v->relBoundRatio, v->pwRelBoundRatio, partial, 1);
				SZ_compress_args_float_ps(&buffer2, (float*)v->data, v->r5, v->r4, v->r3, v->r2, v->r1, &outSize_phase2, v->errBoundMode, v->absErrBound, v->relBoundRatio, v->pwRelBoundRatio, partial, 2);
				
				outSize_[i] = sizeof(size_t) + outSize_phase1 + outSize_phase2;

				unsigned char* p = compressBuffer[i] = (unsigned char*)malloc(outSize_[i]);
				sizeToBytes(p, outSize_phase1);
				p += sizeof(size_t);
				memcpy(p, buffer1, outSize_phase1);
				p += outSize_phase1;
				memcpy(p, buffer2, outSize_phase2);//these 2 step memcpy may be further optimized
				//how to handle compressBuffer[i]?
				free(buffer1); free(buffer2);
				multisteps->compressionType = 2; //this type means the data is compressed both time-based and snapshot-based.
											  //the former part is done by time and the later is done by snapshot.
				//printf("The current varialbe of the current step %d is compressed by method: %d\n", sz_tsc->currentStep, multisteps->compressionType);

			}
			else{
				SZ_compress_args_float(&(compressBuffer[i]), (float*)v->data, v->r5, v->r4, v->r3, v->r2, v->r1, &outSize_[i], v->errBoundMode, v->absErrBound, v->relBoundRatio, v->pwRelBoundRatio);
				//printf("The current varialbe of the current step %d is compressed by method: %d\n", sz_tsc->currentStep, multisteps->compressionType);
			}
			//need to handle the rest part now

		}
		else if(v->dataType==SZ_DOUBLE)
		{
			SZ_compress_args_double(&(compressBuffer[i]), (double*)v->data, v->r5, v->r4, v->r3, v->r2, v->r1, &outSize_[i], v->errBoundMode, v->absErrBound, v->relBoundRatio, v->pwRelBoundRatio);
		}
		sprintf(metadata_str, "%s:%d,%d,%zu", metadata_str, i, multisteps->lastSnapshotStep, outSize_[i]);
		
		totalSize += outSize_[i];
		v->compressType = multisteps->compressionType;
		v = v->next;
	}

#endif
	
	sprintf(metadata_str, "%s\n", metadata_str);
	fputs(metadata_str, sz_tsc->metadata_file);
	free(metadata_str);
	
	//sizeof(int)==current time step; 2*sizeof(char)+sizeof(size_t)=={compressionType + datatype + compression_data_size}; 
	//sizeof(char)==# variables
	//sihuan added: sizeof(size_t)
	*outSize = sizeof(int) + sizeof(unsigned short) + sizeof(size_t) + totalSize + vset->count*(2*sizeof(unsigned char)+sizeof(size_t));
	*newByteData = (unsigned char*)malloc(*outSize); 
	unsigned char* p = *newByteData;

	intToBytes_bigEndian(p, sz_tsc->currentStep);
	p+=4;
	//shortToBytes(p, vset->count);
	shortToBytes(p, vset->count-1);// sihuan updated: subtract the id varialble: so vset->count-1;
	p+=2;
	sizeToBytes(p, sz_tsc->intersect_size);//sihuan added
	p += sizeof(size_t);//sihuan added
	
	//for(i=0;i<vset->count;i++)
	for(i=0;i<vset->count-1;i++)//sihuan updated: subtract the id variable: so vset->count-1;
	{
		//SZ_Variable* v = vset->header->next; //the code is correct? ed?
		if (i == 0) v = vset->header->next;
	
		*p = (unsigned char)v->compressType; //1 byte
		p++;
		*p = (unsigned char)v->dataType; //1 byte
		p++;
		sizeToBytes(p, outSize_[i]); //size_t
		p += sizeof(size_t);
		//sizeToBytes(p, v->r5); //size_t
		//p += sizeof(size_t);
		//sizeToBytes(p, v->r4); //size_t
		//p += sizeof(size_t);
		//sizeToBytes(p, v->r3); //size_t
		//p += sizeof(size_t);
		//sizeToBytes(p, v->r2); //size_t
		//p += sizeof(size_t);
		//sizeToBytes(p, v->r1); //size_t
		//p += sizeof(size_t);								
		memcpy(p, compressBuffer[i], outSize_[i]); //outSize_[i]
		p += outSize_[i];
		v = v->next; //sihuan corrected?
		//free(compressBuffer[i]); //do we need to free it here? why cannot free it?
	}

	//sihuan output size print
	printf("***********************printing results********************\n");
	for (i = 0; i < vset->count - 1; i++)
		printf("variable %d has compressed size: %zu, compression ratio: %.5f\n", i, outSize_[i], (float)(dataLen*sizeof(float))/(float)outSize_[i]);
	printf("***********************ending printing*********************\n");


	sz_tsc->currentStep ++;	
	free(outSize_);
	
	return SZ_SCES;
}



void SZ_decompress_ts(unsigned char *bytes, size_t byteLength)
{
	if(confparams_dec==NULL)
		confparams_dec = (sz_params*)malloc(sizeof(sz_params));
	memset(confparams_dec, 0, sizeof(sz_params));
	confparams_dec->szMode = SZ_TEMPORAL_COMPRESSION;
	confparams_dec->predictionMode = SZ_PREVIOUS_VALUE_ESTIMATE;
	
	if(exe_params==NULL)
		exe_params = (sz_exedata*)malloc(sizeof(sz_exedata));
	memset(exe_params, 0, sizeof(sz_exedata));
	
	int x = 1;
	char *y = (char*)&x;
	if(*y==1)
		sysEndianType = LITTLE_ENDIAN_SYSTEM;
	else //=0
		sysEndianType = BIG_ENDIAN_SYSTEM;
	
	int i = 0;
	size_t r5 = 0, r4 = 0, r3 = 0, r2 = 0, r1 = 0;
	unsigned char* q = bytes;
	sz_tsc->currentStep = bytesToInt_bigEndian(q); 
	//sihuan debug
	printf("currentStep is %d\n", sz_tsc->currentStep);
	q += 4;
	unsigned short nbVars = (unsigned short)bytesToShort(q);
	printf("number of variables is %hu\n", nbVars);//sihuan debug
	q += 2;
	size_t intersection_size = bytesToSize(q);
	q += sizeof(size_t);
	printf("the current step intersection_size is: %zu\n", intersection_size);
	
	if(nbVars != sz_varset->count)
	{
		printf("Error: the number of variables in the compressed data file is inconsistent with the registered # variables.\n");
		printf("Specifically, nbVars = %d, sz_varset->count = %d\n", nbVars, sz_varset->count);
		return;
	}
	
	float *newFloatData = NULL;
	double *newDoubleData = NULL;	
	
	SZ_Variable* p = sz_varset->header->next; // p is pointed to the first variable.
	for(i=0;i<sz_varset->count;i++)
	{
		multisteps = p->multisteps;
		r5 = p->r5;
		r4 = p->r4;
		r3 = p->r3;
		r2 = p->r2;
		r1 = p->r1;
		size_t dataLen = computeDataLength(r5, r4, r3, r2, r1);		
		multisteps->compressionType = *(q++);
		//sihuan debug
		printf("the decompressed compresion type is %d\n", multisteps->compressionType);
		unsigned char dataType = *(q++);
		//sihuan debug
		printf("the decompressed data type is %d\n", SZ_FLOAT);
		size_t cmpSize = bytesToSize(q);
		q += sizeof(size_t);
		unsigned char* cmpBytes = q;
		switch(dataType)
		{
		case SZ_FLOAT:
				if (multisteps->compressionType == 0 || multisteps->compressionType == 1){ //pure snapshot/time compression
					SZ_decompress_args_float(&newFloatData, r5, r4, r3, r2, r1, cmpBytes, cmpSize);
					memcpy(p->data, newFloatData, dataLen*sizeof(float));
					free(newFloatData);
					q += cmpSize;
				}
				else if (multisteps->compressionType == 2){ //mix time, snapshot compression
					size_t tmp_cmpSize1 = bytesToSize(cmpBytes);
					//cmpBytes += sizeof(size_t);
					float* newFloatData1 = NULL;
					float* newFloatData2 = NULL;
					SZ_decompress_args_float_ps(&newFloatData1, r5, r4, r3, r2, r1, cmpBytes+sizeof(size_t), tmp_cmpSize1, 1, intersection_size);
					printf("finish the first decompression phase\n");
					//cmpBytes += tmp_cmpSize;
					size_t tmp_cmpSize2 = cmpSize - tmp_cmpSize1 - sizeof(size_t);
					SZ_decompress_args_float_ps(&newFloatData2, r5, r4, r3, r2, r1, cmpBytes+sizeof(size_t)+tmp_cmpSize1, tmp_cmpSize2, 2, dataLen - intersection_size);
					printf("finish the second decompression phase\n");
					//q += tmp_cmpSize;
					memcpy(p->data, newFloatData1, intersection_size*sizeof(float));
					memcpy(p->data+intersection_size*sizeof(float), newFloatData2, (dataLen-intersection_size)*sizeof(float));
					free(newFloatData1);
					free(newFloatData2);
					q += cmpSize;
				}
				else printf("Only snapshot, time or spatial-time decompression types are supported.\n");
				break;
		case SZ_DOUBLE:
				SZ_decompress_args_double(&newDoubleData, r5, r4, r3, r2, r1, cmpBytes, cmpSize);
				memcpy(p->data, newDoubleData, dataLen*sizeof(double));
				free(newDoubleData);
				break;
		default:
				printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
				return;	
		}
		
		//q += cmpSize;
		p = p->next;
	}
}





void SZ_decompress_ts_vlct(unsigned char *bytes, size_t byteLength)
{
	TheCurVarNum = -1; //each time call this function, this global var should be reset
	if(confparams_dec==NULL)
		confparams_dec = (sz_params*)malloc(sizeof(sz_params));
	memset(confparams_dec, 0, sizeof(sz_params));
	confparams_dec->szMode = SZ_TEMPORAL_COMPRESSION;
	confparams_dec->predictionMode = SZ_PREVIOUS_VALUE_ESTIMATE;
	
	if(exe_params==NULL)
		exe_params = (sz_exedata*)malloc(sizeof(sz_exedata));
	memset(exe_params, 0, sizeof(sz_exedata));
	
	int x = 1;
	char *y = (char*)&x;
	if(*y==1)
		sysEndianType = LITTLE_ENDIAN_SYSTEM;
	else //=0
		sysEndianType = BIG_ENDIAN_SYSTEM;
	
	int i = 0;
	size_t r5 = 0, r4 = 0, r3 = 0, r2 = 0, r1 = 0;
	unsigned char* q = bytes;
	sz_tsc->currentStep = bytesToInt_bigEndian(q); 
	//sihuan debug
	printf("currentStep is %d\n", sz_tsc->currentStep);
	q += 4;
	unsigned short nbVars = (unsigned short)bytesToShort(q);
	printf("number of variables is %hu\n", nbVars);//sihuan debug
	q += 2;
	size_t intersection_size = bytesToSize(q);
	q += sizeof(size_t);
	printf("the current step intersection_size is: %zu\n", intersection_size);
	
	if(nbVars != sz_varset->count)
	{
		printf("Error: the number of variables in the compressed data file is inconsistent with the registered # variables.\n");
		printf("Specifically, nbVars = %d, sz_varset->count = %d\n", nbVars, sz_varset->count);
		return;
	}
	
	float *newFloatData = NULL;
	double *newDoubleData = NULL;	
	
	SZ_Variable* p = sz_varset->header->next; // p is pointed to the first variable.
	for(i=0;i<sz_varset->count;i++)
	{
		TheCurVarNum++;//sihuan added
		v_global = p;
		multisteps = p->multisteps;
		r5 = p->r5;
		r4 = p->r4;
		r3 = p->r3;
		r2 = p->r2;
		r1 = p->r1;
		size_t dataLen = computeDataLength(r5, r4, r3, r2, r1);		
		multisteps->compressionType = *(q++);
		//sihuan debug
		printf("the decompressed compresion type is %d\n", multisteps->compressionType);
		unsigned char dataType = *(q++);
		//sihuan debug
		printf("the decompressed data type is %d\n", SZ_FLOAT);
		size_t cmpSize = bytesToSize(q);
		q += sizeof(size_t);
		unsigned char* cmpBytes = q;
		switch(dataType)
		{
		case SZ_FLOAT:
				if (multisteps->compressionType == 0 || multisteps->compressionType == 1){ //pure snapshot/time compression
					SZ_decompress_args_float(&newFloatData, r5, r4, r3, r2, r1, cmpBytes, cmpSize);
					memcpy(p->data, newFloatData, dataLen*sizeof(float));
					free(newFloatData);
					q += cmpSize;
				}
				else if (multisteps->compressionType == 2){ //mix time, snapshot compression
					size_t tmp_cmpSize1 = bytesToSize(cmpBytes);
					//cmpBytes += sizeof(size_t);
					float* newFloatData1 = NULL;
					float* newFloatData2 = NULL;
					SZ_decompress_args_float_ps(&newFloatData1, r5, r4, r3, r2, r1, cmpBytes+sizeof(size_t), tmp_cmpSize1, 1, intersection_size);
					printf("finish the first decompression phase\n");
					//cmpBytes += tmp_cmpSize;
					size_t tmp_cmpSize2 = cmpSize - tmp_cmpSize1 - sizeof(size_t);
					SZ_decompress_args_float_ps(&newFloatData2, r5, r4, r3, r2, r1, cmpBytes+sizeof(size_t)+tmp_cmpSize1, tmp_cmpSize2, 2, dataLen - intersection_size);
					printf("finish the second decompression phase\n");
					//q += tmp_cmpSize;
					memcpy(p->data, newFloatData1, intersection_size*sizeof(float));
					memcpy(p->data+intersection_size*sizeof(float), newFloatData2, (dataLen-intersection_size)*sizeof(float));
					free(newFloatData1);
					free(newFloatData2);
					q += cmpSize;
				}
				else printf("Only snapshot, time or spatial-time decompression types are supported.\n");
				break;
		case SZ_DOUBLE:
				SZ_decompress_args_double(&newDoubleData, r5, r4, r3, r2, r1, cmpBytes, cmpSize);
				memcpy(p->data, newDoubleData, dataLen*sizeof(double));
				free(newDoubleData);
				break;
		default:
				printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
				return;	
		}
		
		//q += cmpSize;
		p = p->next;
	}
}


#endif


void SZ_Finalize()
{
#ifdef HAVE_TIMECMPR		
	if(sz_varset!=NULL)
		SZ_freeVarSet(SZ_MAINTAIN_VAR_DATA);
#endif

	if(confparams_dec!=NULL)
	{
		free(confparams_dec);
		confparams_dec = NULL;
	}
	if(confparams_cpr!=NULL)
	{
		free(confparams_cpr);
		confparams_cpr = NULL;
	}	
	if(exe_params!=NULL)
	{
		free(exe_params);
		exe_params = NULL;
	}
	
#ifdef HAVE_TIMECMPR	
	if(sz_tsc!=NULL && sz_tsc->metadata_file!=NULL)
		fclose(sz_tsc->metadata_file);
#endif
}
