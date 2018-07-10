/**
 *  @file szd_float.c
 *  @author Sheng Di, Dingwen Tao, Xin Liang
 *  @date Aug, 2018
 *  @brief 
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "szd_float.h"
#include "TightDataPointStorageF.h"
#include "sz.h"
#include "Huffman.h"
#include "szd_float_pwr.h"
#include "szd_float_ts.h"

/**
 * 
 * 
 * @return status SUCCESSFUL (SZ_SCES) or not (other error codes) f
 * */
int SZ_decompress_args_float(float** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, size_t cmpSize)
{
	int status = SZ_SCES;
	size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
	
	//unsigned char* tmpBytes;
	size_t targetUncompressSize = dataLength <<2; //i.e., *4
	//tmpSize must be "much" smaller than dataLength
	size_t i, tmpSize = 8+MetaDataByteLength+exe_params->SZ_SIZE_TYPE;
	unsigned char* szTmpBytes;	
	
	if(cmpSize!=8+4+MetaDataByteLength && cmpSize!=8+8+MetaDataByteLength) //4,8 means two posibilities of SZ_SIZE_TYPE
	{
		int isZlib = isZlibFormat(cmpBytes[0], cmpBytes[1]);
		if(confparams_dec->szMode!=SZ_TEMPORAL_COMPRESSION)
		{
			if(isZlib)
				confparams_dec->szMode = SZ_BEST_COMPRESSION;
			else
				confparams_dec->szMode = SZ_BEST_SPEED;			
		}
		
		if(confparams_dec->szMode==SZ_BEST_SPEED)
		{
			tmpSize = cmpSize;
			szTmpBytes = cmpBytes;	
		}
		else if(confparams_dec->szMode==SZ_BEST_COMPRESSION || confparams_dec->szMode==SZ_DEFAULT_COMPRESSION || confparams_dec->szMode==SZ_TEMPORAL_COMPRESSION)
		{
			if(targetUncompressSize<MIN_ZLIB_DEC_ALLOMEM_BYTES) //Considering the minimum size
				targetUncompressSize = MIN_ZLIB_DEC_ALLOMEM_BYTES; 
			tmpSize = zlib_uncompress5(cmpBytes, (unsigned long)cmpSize, &szTmpBytes, (unsigned long)targetUncompressSize+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE);//		(unsigned long)targetUncompressSize+8: consider the total length under lossless compression mode is actually 3+4+1+targetUncompressSize
			//szTmpBytes = (unsigned char*)malloc(sizeof(unsigned char)*tmpSize);
			//memcpy(szTmpBytes, tmpBytes, tmpSize);
			//free(tmpBytes); //release useless memory		
		}
		else
		{
			printf("Wrong value of confparams_dec->szMode in the double compressed bytes.\n");
			status = SZ_MERR;
			return status;
		}	
	}
	else
		szTmpBytes = cmpBytes;
	//TODO: convert szTmpBytes to data array.
	TightDataPointStorageF* tdps;
	int errBoundMode = new_TightDataPointStorageF_fromFlatBytes(&tdps, szTmpBytes, tmpSize);
	
	//writeByteData(tdps->typeArray, tdps->typeArray_size, "decompress-typebytes.tbt");
	int dim = computeDimension(r5,r4,r3,r2,r1);	
	int floatSize = sizeof(float);
	if(tdps->isLossless)
	{
		*newData = (float*)malloc(floatSize*dataLength);
		if(sysEndianType==BIG_ENDIAN_SYSTEM)
		{
			memcpy(*newData, szTmpBytes+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE, dataLength*floatSize);
		}
		else
		{
			unsigned char* p = szTmpBytes+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE;
			for(i=0;i<dataLength;i++,p+=floatSize)
				(*newData)[i] = bytesToFloat(p);
		}		
	}
	else 
	{
		if(tdps->raBytes_size > 0) //v2.0
		{
			if (dim == 1)
				getSnapshotData_float_1D(newData,r1,tdps, errBoundMode);
			else if(dim == 2)
				decompressDataSeries_float_2D_nonblocked_with_blocked_regression(newData, r2, r1, tdps->raBytes);
			else if(dim == 3)
				decompressDataSeries_float_3D_nonblocked_with_blocked_regression(newData, r3, r2, r1, tdps->raBytes);
			else if(dim == 4)
				decompressDataSeries_float_3D_nonblocked_with_blocked_regression(newData, r4*r3, r2, r1, tdps->raBytes);
			else
			{
				printf("Error: currently support only at most 4 dimensions!\n");
				status = SZ_DERR;
			}	
		}
		else //1.4.13
		{
			if (dim == 1)
				getSnapshotData_float_1D(newData,r1,tdps, errBoundMode);
			else if (dim == 2)
				getSnapshotData_float_2D(newData,r2,r1,tdps, errBoundMode);
			else if (dim == 3)
				getSnapshotData_float_3D(newData,r3,r2,r1,tdps, errBoundMode);
			else if (dim == 4)
				getSnapshotData_float_4D(newData,r4,r3,r2,r1,tdps, errBoundMode);
			else
			{
				printf("Error: currently support only at most 4 dimensions!\n");
				status = SZ_DERR;
			}			
		}
	}
	free_TightDataPointStorageF2(tdps);
	if(confparams_dec->szMode!=SZ_BEST_SPEED && cmpSize!=8+MetaDataByteLength+exe_params->SZ_SIZE_TYPE)
		free(szTmpBytes);
	return status;
}

void decompressDataSeries_float_1D(float** data, size_t dataSeriesLength, TightDataPointStorageF* tdps) 
{
	updateQuantizationInfo(tdps->intervals);
	size_t i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
								// in resiMidBits, p is to track the
								// byte_index of resiMidBits, l is for
								// leadNum
	unsigned char* leadNum;
	double interval = tdps->realPrecision*2;
	
	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (float*)malloc(sizeof(float)*dataSeriesLength);

	int* type = (int*)malloc(dataSeriesLength*sizeof(int));
	
	HuffmanTree* huffmanTree = createHuffmanTree(tdps->stateNum);
	decode_withTree(huffmanTree, tdps->typeArray, dataSeriesLength, type);
	SZ_ReleaseHuffman(huffmanTree);	

	unsigned char preBytes[4];
	unsigned char curBytes[4];
	
	memset(preBytes, 0, 4);

	size_t curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits; 
	unsigned char leadingNum;	
	float medianValue, exactData, predValue;
	
	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;
	
	int type_;
	for (i = 0; i < dataSeriesLength; i++) {	
		type_ = type[i];
		switch (type_) {
		case 0:
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data	
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}
			
			exactData = bytesToFloat(curBytes);
			(*data)[i] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
			break;
		default:
			//predValue = 2 * (*data)[i-1] - (*data)[i-2];
			predValue = (*data)[i-1];
			(*data)[i] = predValue + (type_-exe_params->intvRadius)*interval;
			break;
		}
		//printf("%.30G\n",(*data)[i]);
	}
	
#ifdef HAVE_TIMECMPR	
	if(confparams_dec->szMode == SZ_TEMPORAL_COMPRESSION)
		memcpy(multisteps->hist_data, (*data), dataSeriesLength*sizeof(float));
#endif	
	
	free(leadNum);
	free(type);
	return;
}

void decompressDataSeries_float_2D(float** data, size_t r1, size_t r2, TightDataPointStorageF* tdps) 
{
	updateQuantizationInfo(tdps->intervals);
	//printf("tdps->intervals=%d, exe_params->intvRadius=%d\n", tdps->intervals, exe_params->intvRadius);
	
	size_t j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
	// in resiMidBits, p is to track the
	// byte_index of resiMidBits, l is for
	// leadNum
	size_t dataSeriesLength = r1*r2;
	//	printf ("%d %d\n", r1, r2);

	unsigned char* leadNum;
	double realPrecision = tdps->realPrecision;

	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (float*)malloc(sizeof(float)*dataSeriesLength);

	int* type = (int*)malloc(dataSeriesLength*sizeof(int));

	HuffmanTree* huffmanTree = createHuffmanTree(tdps->stateNum);
	decode_withTree(huffmanTree, tdps->typeArray, dataSeriesLength, type);
	SZ_ReleaseHuffman(huffmanTree);	

	unsigned char preBytes[4];
	unsigned char curBytes[4];

	memset(preBytes, 0, 4);

	size_t curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits; 
	unsigned char leadingNum;	
	float medianValue, exactData;
	int type_;

	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;
	
	float pred1D, pred2D;
	size_t ii, jj;

	/* Process Row-0, data 0 */

	// compute resiBits
	resiBits = 0;
	if (resiBitsLength != 0) {
		int kMod8 = k % 8;
		int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
		if (rightMovSteps > 0) {
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
		} else if (rightMovSteps < 0) {
			int code1 = getLeftMovingCode(kMod8);
			int code2 = getRightMovingCode(kMod8, resiBitsLength);
			int leftMovSteps = -rightMovSteps;
			rightMovSteps = 8 - leftMovSteps;
			resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
			p++;
			resiBits = resiBits
					| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
		} else // rightMovSteps == 0
		{
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code);
			p++;
		}
		k += resiBitsLength;
	}

	// recover the exact data
	memset(curBytes, 0, 4);
	leadingNum = leadNum[l++];
	memcpy(curBytes, preBytes, leadingNum);
	for (j = leadingNum; j < reqBytesLength; j++)
		curBytes[j] = tdps->exactMidBytes[curByteIndex++];
	if (resiBitsLength != 0) {
		unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
		curBytes[reqBytesLength] = resiByte;
	}

	exactData = bytesToFloat(curBytes);
	(*data)[0] = exactData + medianValue;
	memcpy(preBytes,curBytes,4);

	/* Process Row-0, data 1 */
	type_ = type[1]; 
	if (type_ != 0)
	{
		pred1D = (*data)[0];
		(*data)[1] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
	}
	else
	{
		// compute resiBits
		resiBits = 0;
		if (resiBitsLength != 0) {
			int kMod8 = k % 8;
			int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
			if (rightMovSteps > 0) {
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
			} else if (rightMovSteps < 0) {
				int code1 = getLeftMovingCode(kMod8);
				int code2 = getRightMovingCode(kMod8, resiBitsLength);
				int leftMovSteps = -rightMovSteps;
				rightMovSteps = 8 - leftMovSteps;
				resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
				p++;
				resiBits = resiBits
						| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
			} else // rightMovSteps == 0
			{
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code);
				p++;
			}
			k += resiBitsLength;
		}

		// recover the exact data
		memset(curBytes, 0, 4);
		leadingNum = leadNum[l++];
		memcpy(curBytes, preBytes, leadingNum);
		for (j = leadingNum; j < reqBytesLength; j++)
			curBytes[j] = tdps->exactMidBytes[curByteIndex++];
		if (resiBitsLength != 0) {
			unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
			curBytes[reqBytesLength] = resiByte;
		}

		exactData = bytesToFloat(curBytes);
		(*data)[1] = exactData + medianValue;
		memcpy(preBytes,curBytes,4);
	}

	/* Process Row-0, data 2 --> data r2-1 */
	for (jj = 2; jj < r2; jj++)
	{
		type_ = type[jj];
		if (type_ != 0)
		{
			pred1D = 2*(*data)[jj-1] - (*data)[jj-2];				
			(*data)[jj] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[jj] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}
	}

	size_t index;
	/* Process Row-1 --> Row-r1-1 */
	for (ii = 1; ii < r1; ii++)
	{
		/* Process row-ii data 0 */
		index = ii*r2;

		type_ = type[index];
		if (type_ != 0)
		{
			pred1D = (*data)[index-r2];		
			(*data)[index] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}

		/* Process row-ii data 1 --> r2-1*/
		for (jj = 1; jj < r2; jj++)
		{
			index = ii*r2+jj;
			pred2D = (*data)[index-1] + (*data)[index-r2] - (*data)[index-r2-1];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}
		}
	}

#ifdef HAVE_TIMECMPR	
	if(confparams_dec->szMode == SZ_TEMPORAL_COMPRESSION)
		memcpy(multisteps->hist_data, (*data), dataSeriesLength*sizeof(float));
#endif	

	free(leadNum);
	free(type);
	return;
}

void decompressDataSeries_float_3D(float** data, size_t r1, size_t r2, size_t r3, TightDataPointStorageF* tdps) 
{
	updateQuantizationInfo(tdps->intervals);
	size_t j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
	// in resiMidBits, p is to track the
	// byte_index of resiMidBits, l is for
	// leadNum
	size_t dataSeriesLength = r1*r2*r3;
	size_t r23 = r2*r3;
	unsigned char* leadNum;
	double realPrecision = tdps->realPrecision;

	//TODO
	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (float*)malloc(sizeof(float)*dataSeriesLength);
	int* type = (int*)malloc(dataSeriesLength*sizeof(int));

	HuffmanTree* huffmanTree = createHuffmanTree(tdps->stateNum);
	decode_withTree(huffmanTree, tdps->typeArray, dataSeriesLength, type);
	SZ_ReleaseHuffman(huffmanTree);	

	unsigned char preBytes[4];
	unsigned char curBytes[4];

	memset(preBytes, 0, 4);
	size_t curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits;
	unsigned char leadingNum;
	float medianValue, exactData;
	int type_;

	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;
	
	float pred1D, pred2D, pred3D;
	size_t ii, jj, kk;

	///////////////////////////	Process layer-0 ///////////////////////////
	/* Process Row-0 data 0*/
	// compute resiBits
	resiBits = 0;
	if (resiBitsLength != 0) {
		int kMod8 = k % 8;
		int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
		if (rightMovSteps > 0) {
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
		} else if (rightMovSteps < 0) {
			int code1 = getLeftMovingCode(kMod8);
			int code2 = getRightMovingCode(kMod8, resiBitsLength);
			int leftMovSteps = -rightMovSteps;
			rightMovSteps = 8 - leftMovSteps;
			resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
			p++;
			resiBits = resiBits
					| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
		} else // rightMovSteps == 0
		{
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code);
			p++;
		}
		k += resiBitsLength;
	}

	// recover the exact data
	memset(curBytes, 0, 4);
	leadingNum = leadNum[l++];
	memcpy(curBytes, preBytes, leadingNum);
	for (j = leadingNum; j < reqBytesLength; j++)
		curBytes[j] = tdps->exactMidBytes[curByteIndex++];
	if (resiBitsLength != 0) {
		unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
		curBytes[reqBytesLength] = resiByte;
	}
	exactData = bytesToFloat(curBytes);
	(*data)[0] = exactData + medianValue;
	memcpy(preBytes,curBytes,4);

	/* Process Row-0, data 1 */
	pred1D = (*data)[0];

	type_ = type[1];
	if (type_ != 0)
	{
		(*data)[1] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
	}
	else
	{
		// compute resiBits
		resiBits = 0;
		if (resiBitsLength != 0) {
			int kMod8 = k % 8;
			int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
			if (rightMovSteps > 0) {
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
			} else if (rightMovSteps < 0) {
				int code1 = getLeftMovingCode(kMod8);
				int code2 = getRightMovingCode(kMod8, resiBitsLength);
				int leftMovSteps = -rightMovSteps;
				rightMovSteps = 8 - leftMovSteps;
				resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
				p++;
				resiBits = resiBits
						| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
			} else // rightMovSteps == 0
			{
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code);
				p++;
			}
			k += resiBitsLength;
		}

		// recover the exact data
		memset(curBytes, 0, 4);
		leadingNum = leadNum[l++];
		memcpy(curBytes, preBytes, leadingNum);
		for (j = leadingNum; j < reqBytesLength; j++)
			curBytes[j] = tdps->exactMidBytes[curByteIndex++];
		if (resiBitsLength != 0) {
			unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
			curBytes[reqBytesLength] = resiByte;
		}

		exactData = bytesToFloat(curBytes);
		(*data)[1] = exactData + medianValue;
		memcpy(preBytes,curBytes,4);
	}
	/* Process Row-0, data 2 --> data r3-1 */
	for (jj = 2; jj < r3; jj++)
	{
		pred1D = 2*(*data)[jj-1] - (*data)[jj-2];

		type_ = type[jj];
		if (type_ != 0)
		{
			(*data)[jj] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[jj] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}
	}

	size_t index;
	/* Process Row-1 --> Row-r2-1 */
	for (ii = 1; ii < r2; ii++)
	{
		/* Process row-ii data 0 */
		index = ii*r3;
		pred1D = (*data)[index-r3];

		type_ = type[index];
		if (type_ != 0)
		{
			(*data)[index] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}

		/* Process row-ii data 1 --> r3-1*/
		for (jj = 1; jj < r3; jj++)
		{
			index = ii*r3+jj;
			pred2D = (*data)[index-1] + (*data)[index-r3] - (*data)[index-r3-1];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}
		}
	}

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (kk = 1; kk < r1; kk++)
	{
		/* Process Row-0 data 0*/
		index = kk*r23;
		pred1D = (*data)[index-r23];

		type_ = type[index];
		if (type_ != 0)
		{
			(*data)[index] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}

		/* Process Row-0 data 1 --> data r3-1 */
		for (jj = 1; jj < r3; jj++)
		{
			index = kk*r23+jj;
			pred2D = (*data)[index-1] + (*data)[index-r23] - (*data)[index-r23-1];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}
		}

		/* Process Row-1 --> Row-r2-1 */
		for (ii = 1; ii < r2; ii++)
		{
			/* Process Row-i data 0 */
			index = kk*r23 + ii*r3;
			pred2D = (*data)[index-r3] + (*data)[index-r23] - (*data)[index-r23-r3];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (jj = 1; jj < r3; jj++)
			{
				index = kk*r23 + ii*r3 + jj;
				pred3D = (*data)[index-1] + (*data)[index-r3] + (*data)[index-r23]
					- (*data)[index-r3-1] - (*data)[index-r23-r3] - (*data)[index-r23-1] + (*data)[index-r23-r3-1];

				type_ = type[index];
				if (type_ != 0)
				{
					(*data)[index] = pred3D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
				}
				else
				{
					// compute resiBits
					resiBits = 0;
					if (resiBitsLength != 0) {
						int kMod8 = k % 8;
						int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
						if (rightMovSteps > 0) {
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
						} else if (rightMovSteps < 0) {
							int code1 = getLeftMovingCode(kMod8);
							int code2 = getRightMovingCode(kMod8, resiBitsLength);
							int leftMovSteps = -rightMovSteps;
							rightMovSteps = 8 - leftMovSteps;
							resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
							p++;
							resiBits = resiBits
									| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
						} else // rightMovSteps == 0
						{
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code);
							p++;
						}
						k += resiBitsLength;
					}

					// recover the exact data
					memset(curBytes, 0, 4);
					leadingNum = leadNum[l++];
					memcpy(curBytes, preBytes, leadingNum);
					for (j = leadingNum; j < reqBytesLength; j++)
						curBytes[j] = tdps->exactMidBytes[curByteIndex++];
					if (resiBitsLength != 0) {
						unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
						curBytes[reqBytesLength] = resiByte;
					}

					exactData = bytesToFloat(curBytes);
					(*data)[index] = exactData + medianValue;
					memcpy(preBytes,curBytes,4);
				}
			}
		}
	}
	
#ifdef HAVE_TIMECMPR	
	if(confparams_dec->szMode == SZ_TEMPORAL_COMPRESSION)
		memcpy(multisteps->hist_data, (*data), dataSeriesLength*sizeof(float));
#endif		

	free(leadNum);
	free(type);
	return;
}


void decompressDataSeries_float_4D(float** data, size_t r1, size_t r2, size_t r3, size_t r4, TightDataPointStorageF* tdps)
{
	updateQuantizationInfo(tdps->intervals);
	size_t j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
	// in resiMidBits, p is to track the
	// byte_index of resiMidBits, l is for
	// leadNum
	size_t dataSeriesLength = r1*r2*r3*r4;
	size_t r234 = r2*r3*r4;
	size_t r34 = r3*r4;
//	printf ("%d %d %d %d\n", r1, r2, r3, r4);
	unsigned char* leadNum;
	double realPrecision = tdps->realPrecision;

	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (float*)malloc(sizeof(float)*dataSeriesLength);
	int* type = (int*)malloc(dataSeriesLength*sizeof(int));

	HuffmanTree* huffmanTree = createHuffmanTree(tdps->stateNum);
	decode_withTree(huffmanTree, tdps->typeArray, dataSeriesLength, type);
	SZ_ReleaseHuffman(huffmanTree);	

	unsigned char preBytes[4];
	unsigned char curBytes[4];

	memset(preBytes, 0, 4);
	size_t curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits;
	unsigned char leadingNum;
	float medianValue, exactData;
	int type_;

	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;

	float pred1D, pred2D, pred3D;
	size_t ii, jj, kk, ll;
	size_t index;

	for (ll = 0; ll < r1; ll++)
	{

		///////////////////////////	Process layer-0 ///////////////////////////
		/* Process Row-0 data 0*/
		index = ll*r234;

		// compute resiBits
		resiBits = 0;
		if (resiBitsLength != 0) {
			int kMod8 = k % 8;
			int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
			if (rightMovSteps > 0) {
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
			} else if (rightMovSteps < 0) {
				int code1 = getLeftMovingCode(kMod8);
				int code2 = getRightMovingCode(kMod8, resiBitsLength);
				int leftMovSteps = -rightMovSteps;
				rightMovSteps = 8 - leftMovSteps;
				resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
				p++;
				resiBits = resiBits
						| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
			} else // rightMovSteps == 0
			{
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code);
				p++;
			}
			k += resiBitsLength;
		}

		// recover the exact data
		memset(curBytes, 0, 4);
		leadingNum = leadNum[l++];
		memcpy(curBytes, preBytes, leadingNum);
		for (j = leadingNum; j < reqBytesLength; j++)
			curBytes[j] = tdps->exactMidBytes[curByteIndex++];
		if (resiBitsLength != 0) {
			unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
			curBytes[reqBytesLength] = resiByte;
		}
		exactData = bytesToFloat(curBytes);
		(*data)[index] = exactData + medianValue;
		memcpy(preBytes,curBytes,4);

		/* Process Row-0, data 1 */
		index = ll*r234+1;

		pred1D = (*data)[index-1];

		type_ = type[index];
		if (type_ != 0)
		{
			(*data)[index] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}

		/* Process Row-0, data 2 --> data r4-1 */
		for (jj = 2; jj < r4; jj++)
		{
			index = ll*r234+jj;

			pred1D = 2*(*data)[index-1] - (*data)[index-2];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}
		}

		/* Process Row-1 --> Row-r3-1 */
		for (ii = 1; ii < r3; ii++)
		{
			/* Process row-ii data 0 */
			index = ll*r234+ii*r4;

			pred1D = (*data)[index-r4];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}

			/* Process row-ii data 1 --> r4-1*/
			for (jj = 1; jj < r4; jj++)
			{
				index = ll*r234+ii*r4+jj;

				pred2D = (*data)[index-1] + (*data)[index-r4] - (*data)[index-r4-1];

				type_ = type[index];
				if (type_ != 0)
				{
					(*data)[index] = pred2D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
				}
				else
				{
					// compute resiBits
					resiBits = 0;
					if (resiBitsLength != 0) {
						int kMod8 = k % 8;
						int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
						if (rightMovSteps > 0) {
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
						} else if (rightMovSteps < 0) {
							int code1 = getLeftMovingCode(kMod8);
							int code2 = getRightMovingCode(kMod8, resiBitsLength);
							int leftMovSteps = -rightMovSteps;
							rightMovSteps = 8 - leftMovSteps;
							resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
							p++;
							resiBits = resiBits
									| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
						} else // rightMovSteps == 0
						{
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code);
							p++;
						}
						k += resiBitsLength;
					}

					// recover the exact data
					memset(curBytes, 0, 4);
					leadingNum = leadNum[l++];
					memcpy(curBytes, preBytes, leadingNum);
					for (j = leadingNum; j < reqBytesLength; j++)
						curBytes[j] = tdps->exactMidBytes[curByteIndex++];
					if (resiBitsLength != 0) {
						unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
						curBytes[reqBytesLength] = resiByte;
					}

					exactData = bytesToFloat(curBytes);
					(*data)[index] = exactData + medianValue;
					memcpy(preBytes,curBytes,4);
				}
			}
		}

		///////////////////////////	Process layer-1 --> layer-r2-1 ///////////////////////////

		for (kk = 1; kk < r2; kk++)
		{
			/* Process Row-0 data 0*/
			index = ll*r234+kk*r34;

			pred1D = (*data)[index-r34];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}

			/* Process Row-0 data 1 --> data r4-1 */
			for (jj = 1; jj < r4; jj++)
			{
				index = ll*r234+kk*r34+jj;

				pred2D = (*data)[index-1] + (*data)[index-r34] - (*data)[index-r34-1];

				type_ = type[index];
				if (type_ != 0)
				{
					(*data)[index] = pred2D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
				}
				else
				{
					// compute resiBits
					resiBits = 0;
					if (resiBitsLength != 0) {
						int kMod8 = k % 8;
						int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
						if (rightMovSteps > 0) {
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
						} else if (rightMovSteps < 0) {
							int code1 = getLeftMovingCode(kMod8);
							int code2 = getRightMovingCode(kMod8, resiBitsLength);
							int leftMovSteps = -rightMovSteps;
							rightMovSteps = 8 - leftMovSteps;
							resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
							p++;
							resiBits = resiBits
									| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
						} else // rightMovSteps == 0
						{
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code);
							p++;
						}
						k += resiBitsLength;
					}

					// recover the exact data
					memset(curBytes, 0, 4);
					leadingNum = leadNum[l++];
					memcpy(curBytes, preBytes, leadingNum);
					for (j = leadingNum; j < reqBytesLength; j++)
						curBytes[j] = tdps->exactMidBytes[curByteIndex++];
					if (resiBitsLength != 0) {
						unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
						curBytes[reqBytesLength] = resiByte;
					}

					exactData = bytesToFloat(curBytes);
					(*data)[index] = exactData + medianValue;
					memcpy(preBytes,curBytes,4);
				}
			}

			/* Process Row-1 --> Row-r3-1 */
			for (ii = 1; ii < r3; ii++)
			{
				/* Process Row-i data 0 */
				index = ll*r234+kk*r34+ii*r4;

				pred2D = (*data)[index-r4] + (*data)[index-r34] - (*data)[index-r34-r4];

				type_ = type[index];
				if (type_ != 0)
				{
					(*data)[index] = pred2D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
				}
				else
				{
					// compute resiBits
					resiBits = 0;
					if (resiBitsLength != 0) {
						int kMod8 = k % 8;
						int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
						if (rightMovSteps > 0) {
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
						} else if (rightMovSteps < 0) {
							int code1 = getLeftMovingCode(kMod8);
							int code2 = getRightMovingCode(kMod8, resiBitsLength);
							int leftMovSteps = -rightMovSteps;
							rightMovSteps = 8 - leftMovSteps;
							resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
							p++;
							resiBits = resiBits
									| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
						} else // rightMovSteps == 0
						{
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code);
							p++;
						}
						k += resiBitsLength;
					}

					// recover the exact data
					memset(curBytes, 0, 4);
					leadingNum = leadNum[l++];
					memcpy(curBytes, preBytes, leadingNum);
					for (j = leadingNum; j < reqBytesLength; j++)
						curBytes[j] = tdps->exactMidBytes[curByteIndex++];
					if (resiBitsLength != 0) {
						unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
						curBytes[reqBytesLength] = resiByte;
					}

					exactData = bytesToFloat(curBytes);
					(*data)[index] = exactData + medianValue;
					memcpy(preBytes,curBytes,4);
				}

				/* Process Row-i data 1 --> data r4-1 */
				for (jj = 1; jj < r4; jj++)
				{
					index = ll*r234+kk*r34+ii*r4+jj;

					pred3D = (*data)[index-1] + (*data)[index-r4] + (*data)[index-r34]
							- (*data)[index-r4-1] - (*data)[index-r34-r4] - (*data)[index-r34-1] + (*data)[index-r34-r4-1];


					type_ = type[index];
					if (type_ != 0)
					{
						(*data)[index] = pred3D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
					}
					else
					{
						// compute resiBits
						resiBits = 0;
						if (resiBitsLength != 0) {
							int kMod8 = k % 8;
							int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
							if (rightMovSteps > 0) {
								int code = getRightMovingCode(kMod8, resiBitsLength);
								resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
							} else if (rightMovSteps < 0) {
								int code1 = getLeftMovingCode(kMod8);
								int code2 = getRightMovingCode(kMod8, resiBitsLength);
								int leftMovSteps = -rightMovSteps;
								rightMovSteps = 8 - leftMovSteps;
								resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
								p++;
								resiBits = resiBits
										| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
							} else // rightMovSteps == 0
							{
								int code = getRightMovingCode(kMod8, resiBitsLength);
								resiBits = (tdps->residualMidBits[p] & code);
								p++;
							}
							k += resiBitsLength;
						}

						// recover the exact data
						memset(curBytes, 0, 4);
						leadingNum = leadNum[l++];
						memcpy(curBytes, preBytes, leadingNum);
						for (j = leadingNum; j < reqBytesLength; j++)
							curBytes[j] = tdps->exactMidBytes[curByteIndex++];
						if (resiBitsLength != 0) {
							unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
							curBytes[reqBytesLength] = resiByte;
						}

						exactData = bytesToFloat(curBytes);
						(*data)[index] = exactData + medianValue;
						memcpy(preBytes,curBytes,4);
					}
				}
			}

		}
	}

//I didn't implement time-based compression for 4D actually. 
//#ifdef HAVE_TIMECMPR	
//	if(confparams_dec->szMode == SZ_TEMPORAL_COMPRESSION)
//		memcpy(multisteps->hist_data, (*data), dataSeriesLength*sizeof(float));
//#endif	

	free(leadNum);
	free(type);
	return;
}

void getSnapshotData_float_1D(float** data, size_t dataSeriesLength, TightDataPointStorageF* tdps, int errBoundMode)
{	
	size_t i;

	if (tdps->allSameData) {
		float value = bytesToFloat(tdps->exactMidBytes);
		*data = (float*)malloc(sizeof(float)*dataSeriesLength);
		for (i = 0; i < dataSeriesLength; i++)
			(*data)[i] = value;
	} else {
		if (tdps->rtypeArray == NULL) {
			if(errBoundMode < PW_REL)
			{
#ifdef HAVE_TIMECMPR				
				if(confparams_dec->szMode == SZ_TEMPORAL_COMPRESSION)
				{
					if(multisteps->compressionType == 0) //snapshot
						decompressDataSeries_float_1D(data, dataSeriesLength, tdps);
					else
						decompressDataSeries_float_1D_ts(data, dataSeriesLength, multisteps, tdps);					
				}
				else
#endif				
					decompressDataSeries_float_1D(data, dataSeriesLength, tdps);
			}
			else 
			{
				//decompressDataSeries_float_1D_pwr(data, dataSeriesLength, tdps);
				decompressDataSeries_float_1D_pwrgroup(data, dataSeriesLength, tdps);
			}
			return;
		} else {
			*data = (float*)malloc(sizeof(float)*dataSeriesLength);
			// insert the reserved values
			//int[] rtypes = TypeManager.convertByteArray2IntArray_fast_1b(
			//		dataSeriesLength, rtypeArray);
			int* rtypes;
			int validLength = computeBitNumRequired(dataSeriesLength);
			decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
			size_t count = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 1)
					(*data)[i] = tdps->reservedValue;
				else
					count++;
			}
			// get the decompressed data
			float* decmpData;
			if(errBoundMode < PW_REL)
				decompressDataSeries_float_1D(&decmpData, dataSeriesLength, tdps);
			else 
				decompressDataSeries_float_1D_pwr(&decmpData, dataSeriesLength, tdps);
			// insert the decompressed data
			size_t k = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 0) {
					(*data)[i] = decmpData[k++];
				}
			}
			free(decmpData);
			free(rtypes);
		}
	}
}

void getSnapshotData_float_2D(float** data, size_t r1, size_t r2, TightDataPointStorageF* tdps, int errBoundMode) 
{
	size_t i;
	size_t dataSeriesLength = r1*r2;
	if (tdps->allSameData) {
		float value = bytesToFloat(tdps->exactMidBytes);
		*data = (float*)malloc(sizeof(float)*dataSeriesLength);
		for (i = 0; i < dataSeriesLength; i++)
			(*data)[i] = value;
	} else {
		if (tdps->rtypeArray == NULL) {
			if(errBoundMode < PW_REL)
			{
#ifdef HAVE_TIMECMPR					
				if(confparams_dec->szMode == SZ_TEMPORAL_COMPRESSION)
				{
					if(multisteps->compressionType == 0)
						decompressDataSeries_float_2D(data, r1, r2, tdps);
					else
						decompressDataSeries_float_1D_ts(data, r1*r2, multisteps, tdps);					
				}
				else
#endif
					decompressDataSeries_float_2D(data, r1, r2, tdps);
			}
			else 
			{
				decompressDataSeries_float_2D_pwr(data, r1, r2, tdps);
			}			

			return;
		} else {
			*data = (float*)malloc(sizeof(float)*dataSeriesLength);
			// insert the reserved values
			//int[] rtypes = TypeManager.convertByteArray2IntArray_fast_1b(
			//		dataSeriesLength, rtypeArray);
			int* rtypes;
			int validLength = computeBitNumRequired(dataSeriesLength);
			decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
			size_t count = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 1)
					(*data)[i] = tdps->reservedValue;
				else
					count++;
			}
			// get the decompressed data
			float* decmpData;
			if(errBoundMode < PW_REL)
				decompressDataSeries_float_2D(&decmpData, r1, r2, tdps);
			else 
				decompressDataSeries_float_2D_pwr(&decmpData, r1, r2, tdps);
			// insert the decompressed data
			size_t k = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 0) {
					(*data)[i] = decmpData[k++];
				}
			}
			free(decmpData);
			free(rtypes);
		}
	}
}

void getSnapshotData_float_3D(float** data, size_t r1, size_t r2, size_t r3, TightDataPointStorageF* tdps, int errBoundMode)
{
	size_t i;
	size_t dataSeriesLength = r1*r2*r3;
	if (tdps->allSameData) {
		float value = bytesToFloat(tdps->exactMidBytes);
		*data = (float*)malloc(sizeof(float)*dataSeriesLength);
		for (i = 0; i < dataSeriesLength; i++)
			(*data)[i] = value;
	} else {
		if (tdps->rtypeArray == NULL) {
			if(errBoundMode < PW_REL)
			{
#ifdef HAVE_TIMECMPR					
				if(confparams_dec->szMode == SZ_TEMPORAL_COMPRESSION)
				{
					if(multisteps->compressionType == 0)
						decompressDataSeries_float_3D(data, r1, r2, r3, tdps);
					else
						decompressDataSeries_float_1D_ts(data, r1*r2*r3, multisteps, tdps);					
				}
				else
#endif				
					decompressDataSeries_float_3D(data, r1, r2, r3, tdps);
			}
			else 
			{
				decompressDataSeries_float_3D_pwr(data, r1, r2, r3, tdps);
			}					
			
			return;
		} else {
			*data = (float*)malloc(sizeof(float)*dataSeriesLength);
			// insert the reserved values
			//int[] rtypes = TypeManager.convertByteArray2IntArray_fast_1b(
			//		dataSeriesLength, rtypeArray);
			int* rtypes;
			int validLength = computeBitNumRequired(dataSeriesLength);
			decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
			size_t count = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 1)
					(*data)[i] = tdps->reservedValue;
				else
					count++;
			}
			// get the decompressed data
			float* decmpData;
			if(errBoundMode < PW_REL)
				decompressDataSeries_float_3D(&decmpData, r1, r2, r3, tdps);
			else 
				decompressDataSeries_float_3D_pwr(&decmpData, r1, r2, r3, tdps);
			// insert the decompressed data
			size_t k = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 0) {
					(*data)[i] = decmpData[k++];
				}
			}
			free(decmpData);
			free(rtypes);
		}
	}
}

void getSnapshotData_float_4D(float** data, size_t r1, size_t r2, size_t r3, size_t r4, TightDataPointStorageF* tdps, int errBoundMode)
{
	size_t i;
	size_t dataSeriesLength = r1*r2*r3*r4;
	if (tdps->allSameData) {
		float value = bytesToFloat(tdps->exactMidBytes);
		*data = (float*)malloc(sizeof(float)*dataSeriesLength);
		for (i = 0; i < dataSeriesLength; i++)
			(*data)[i] = value;
	} else {
		if (tdps->rtypeArray == NULL) {
			if(errBoundMode < PW_REL)
			{
#ifdef HAVE_TIMECMPR					
				if(confparams_dec->szMode == SZ_TEMPORAL_COMPRESSION)
				{
					if(multisteps->compressionType == 0)
						decompressDataSeries_float_4D(data, r1, r2, r3, r4, tdps);
					else
						decompressDataSeries_float_1D_ts(data, r1*r2*r3*r4, multisteps, tdps);					
				}
				else
#endif				
					decompressDataSeries_float_4D(data, r1, r2, r3, r4, tdps);
			}
			else 
			{
				decompressDataSeries_float_3D_pwr(data, r1*r2, r3, r4, tdps);
				//ToDO
				//decompressDataSeries_float_4D_pwr(data, r1, r2, r3, r4, tdps);
			}					
			return;
		} else {
			*data = (float*)malloc(sizeof(float)*dataSeriesLength);
			int* rtypes;
			int validLength = computeBitNumRequired(dataSeriesLength);
			decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
			size_t count = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 1)
					(*data)[i] = tdps->reservedValue;
				else
					count++;
			}
			// get the decompressed data
			float* decmpData;
			if(errBoundMode < PW_REL)
				decompressDataSeries_float_4D(&decmpData, r1, r2, r3, r4, tdps);
			else
				decompressDataSeries_float_3D_pwr(&decmpData, r1*r2, r3, r4, tdps);
				//ToDO
				//decompressDataSeries_float_4D_pwr(&decompData, r1, r2, r3, r4, tdps);
			// insert the decompressed data
			size_t k = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 0) {
					(*data)[i] = decmpData[k++];
				}
			}
			free(decmpData);
			free(rtypes);
		}
	}
}

size_t decompressDataSeries_float_3D_RA_block(float * data, float mean, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, int * type, float * unpredictable_data){

	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;
	// printf("SZ_compress_float_3D_MDQ_RA_block real dim: %d %d %d\n", real_block_dims[0], real_block_dims[1], real_block_dims[2]);
	// fflush(stdout);

	size_t unpredictable_count = 0;
	size_t r1, r2, r3;
	r1 = block_dim_0;
	r2 = block_dim_1;
	r3 = block_dim_2;

	float * cur_data_pos = data;
	float * last_row_pos;
	float pred1D, pred2D, pred3D;
	size_t i, j, k;
	size_t r23 = r2*r3;
	int type_;
	// Process Row-0 data 0
	pred1D = mean;
	type_ = type[0];
	// printf("Type 0 %d, mean %.4f\n", type_, mean);
	if (type_ != 0){
		cur_data_pos[0] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
	}

	/* Process Row-0 data 1*/
	pred1D = cur_data_pos[0];
	type_ = type[1];
	if (type_ != 0){
		cur_data_pos[1] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[1] = unpredictable_data[unpredictable_count ++];
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++){
		pred1D = 2*cur_data_pos[j-1] - cur_data_pos[j-2];
		type_ = type[j];
		if (type_ != 0){
			cur_data_pos[j] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
		}
	}

	last_row_pos = cur_data_pos;
	cur_data_pos += dim1_offset;
	// printf("SZ_compress_float_3D_MDQ_RA_block row 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;	
		pred1D = last_row_pos[0];
		type_ = type[index];
		if (type_ != 0){
			cur_data_pos[0] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			pred2D = cur_data_pos[j-1] + last_row_pos[j] - last_row_pos[j-1];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[j] = pred2D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			// printf("pred2D %.2f cur_data %.2f last_row_data %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], last_row_pos[j], last_row_pos[j-1], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim1_offset;
	}
	cur_data_pos += dim0_offset - r2 * dim1_offset;

	// printf("SZ_compress_float_3D_MDQ_RA_block layer 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);
	// exit(0);

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		// if(idx == 63 && idy == 63 && idz == 63){
		// 	printf("SZ_compress_float_3D_MDQ_RA_block layer %d done, cur_data_pos: %ld\n", k-1, cur_data_pos - data);
		// 	fflush(stdout);
		// }
		/* Process Row-0 data 0*/
		index = k*r23;
		pred1D = cur_data_pos[- dim0_offset];
		type_ = type[index];
		if (type_ != 0){
			cur_data_pos[0] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;
			index ++;
			pred2D = cur_data_pos[j-1] + cur_data_pos[j - dim0_offset] - cur_data_pos[j - 1 - dim0_offset];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[j] = pred2D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			// printf("pred2D %.2f cur_data %.2f %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], cur_data_pos[j - dim0_offset], cur_data_pos[j - 1 - dim0_offset], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim1_offset;

		// if(idx == 63 && idy == 63 && idz == 63){
		// 	printf("SZ_compress_float_3D_MDQ_RA_block layer row 0 done, cur_data_pos: %ld\n", k-1, cur_data_pos - data);
		// 	fflush(stdout);
		// }

	    /* Process Row-1 --> Row-r2-1 */
		for (i = 1; i < r2; i++)
		{
			// if(idx == 63 && idy == 63 && idz == 63){
			// 	printf("SZ_compress_float_3D_MDQ_RA_block layer row %d done, cur_data_pos: %ld\n", i-1, cur_data_pos - data);
			// 	fflush(stdout);
			// }
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			pred2D = last_row_pos[0] + cur_data_pos[- dim0_offset] - last_row_pos[- dim0_offset];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[0] = pred2D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
//				if(k==63&&i==43&&j==27)
//					printf("i=%d\n", i);
				//index = k*r2*r3 + i*r3 + j;			
				index ++;
				pred3D = cur_data_pos[j-1] + last_row_pos[j]+ cur_data_pos[j - dim0_offset] - last_row_pos[j-1] - last_row_pos[j - dim0_offset] - cur_data_pos[j-1 - dim0_offset] + last_row_pos[j-1 - dim0_offset];
				type_ = type[index];
				if (type_ != 0){
					cur_data_pos[j] = pred3D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
				}
				else{
					cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
				}
			}
			last_row_pos = cur_data_pos;
			cur_data_pos += dim1_offset;
		}
		cur_data_pos += dim0_offset - r2 * dim1_offset;
	}

	return unpredictable_count;
}

size_t decompressDataSeries_float_1D_RA_block(float * data, float mean, size_t dim_0, size_t block_dim_0, double realPrecision, int * type, float * unpredictable_data){

	size_t unpredictable_count = 0;
	
	float * cur_data_pos = data;
	size_t type_index = 0;
	int type_;
	float last_over_thres = mean;
	for(size_t i=0; i<block_dim_0; i++){
		type_ = type[type_index];
		if(type_ == 0){
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
			last_over_thres = cur_data_pos[0];
		}
		else{
			cur_data_pos[0] = last_over_thres + 2 * (type_ - exe_params->intvRadius) * realPrecision;
			last_over_thres = cur_data_pos[0];
		}

		type_index ++;
		cur_data_pos ++;
	}

	return unpredictable_count;
}

size_t decompressDataSeries_float_2D_RA_block(float * data, float mean, size_t dim_0, size_t dim_1, size_t block_dim_0, size_t block_dim_1, double realPrecision, int * type, float * unpredictable_data){

	size_t dim0_offset = dim_1;
	// printf("SZ_compress_float_3D_MDQ_RA_block real dim: %d %d %d\n", real_block_dims[0], real_block_dims[1], real_block_dims[2]);
	// fflush(stdout);

	size_t unpredictable_count = 0;
	size_t r1, r2;
	r1 = block_dim_0;
	r2 = block_dim_1;

	float * cur_data_pos = data;
	float * last_row_pos;
	float pred1D, pred2D;
	size_t i, j;
	int type_;
	// Process Row-0 data 0
	pred1D = mean;
	type_ = type[0];
	// printf("Type 0 %d, mean %.4f\n", type_, mean);
	if (type_ != 0){
		cur_data_pos[0] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
	}

	/* Process Row-0 data 1*/
	pred1D = cur_data_pos[0];
	type_ = type[1];
	if (type_ != 0){
		cur_data_pos[1] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[1] = unpredictable_data[unpredictable_count ++];
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r2; j++){
		pred1D = 2*cur_data_pos[j-1] - cur_data_pos[j-2];
		type_ = type[j];
		if (type_ != 0){
			cur_data_pos[j] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
		}
	}

	last_row_pos = cur_data_pos;
	cur_data_pos += dim0_offset;
	// printf("SZ_compress_float_3D_MDQ_RA_block row 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r1; i++)
	{
		/* Process row-i data 0 */
		index = i*r2;	
		type_ = type[index];
		if (type_ != 0){
			pred1D = last_row_pos[0];
			cur_data_pos[0] = pred1D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r2; j++)
		{
			index = i*r2+j;
			pred2D = cur_data_pos[j-1] + last_row_pos[j] - last_row_pos[j-1];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[j] = pred2D + 2 * (type_ - exe_params->intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			// printf("pred2D %.2f cur_data %.2f last_row_data %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], last_row_pos[j], last_row_pos[j-1], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim0_offset;
	}
	return unpredictable_count;
}

void decompressDataSeries_float_2D_nonblocked_with_blocked_regression(float** data, size_t r1, size_t r2, unsigned char* comp_data){

	size_t dim0_offset = r2;
	size_t num_elements = r1 * r2;

	*data = (float*)malloc(sizeof(float)*num_elements);

	unsigned char * comp_data_pos = comp_data;

	size_t block_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += sizeof(int);
	// calculate block dims
	size_t num_x, num_y;
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r2, num_y, block_size);

	size_t split_index_x, split_index_y;
	size_t early_blockcount_x, early_blockcount_y;
	size_t late_blockcount_x, late_blockcount_y;
	SZ_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	SZ_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);

	size_t num_blocks = num_x * num_y;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += sizeof(double);
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += sizeof(int);

	updateQuantizationInfo(intervals);

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += sizeof(int);

	int stateNum = 2*intervals;
	HuffmanTree* huffmanTree = createHuffmanTree(stateNum);
	
	int nodeCount = bytesToInt_bigEndian(comp_data_pos);
	
	node root = reconstruct_HuffTree_from_bytes_anyStates(huffmanTree,comp_data_pos+4, nodeCount);
	comp_data_pos += sizeof(int) + tree_size;

	float mean;
	unsigned char use_mean;
	memcpy(&use_mean, comp_data_pos, sizeof(unsigned char));
	comp_data_pos += sizeof(unsigned char);
	memcpy(&mean, comp_data_pos, sizeof(float));
	comp_data_pos += sizeof(float);
	size_t reg_count = 0;

	unsigned char * indicator;
	size_t indicator_bitlength = (num_blocks - 1)/8 + 1;
	convertByteArray2IntArray_fast_1b(num_blocks, comp_data_pos, indicator_bitlength, &indicator);
	comp_data_pos += indicator_bitlength;
	for(size_t i=0; i<num_blocks; i++){
		if(!indicator[i]) reg_count ++;
	}
	//printf("reg_count: %ld\n", reg_count);

	float* dec_a = NULL, *dec_b = NULL, *dec_c = NULL;
	if(reg_count > 0){
		float * medians = (float *) comp_data_pos;
		float medianValue_a = *medians;
		float medianValue_b = *(medians + 1);
		float medianValue_c = *(medians + 2);
		comp_data_pos += 3 * sizeof(float);

		int * reqLength = (int *) comp_data_pos;
		int reqLength_a = *reqLength;
		int reqLength_b = *(reqLength + 1);
		int reqLength_c = *(reqLength + 2);
		comp_data_pos += 3 * sizeof(int);
		
		//reconstruct leading_number array in the form of meaningful integers....
		size_t leadNumArray_size = (reg_count-1)/4+1;
		
		unsigned char* leadNum_a = NULL, *leadNum_b = NULL, *leadNum_c = NULL;
		
		unsigned char* leadNumArray = comp_data_pos;
		convertByteArray2IntArray_fast_2b(reg_count, leadNumArray, leadNumArray_size, &leadNum_a);
		comp_data_pos += leadNumArray_size;
		
		leadNumArray = comp_data_pos;
		convertByteArray2IntArray_fast_2b(reg_count, leadNumArray, leadNumArray_size, &leadNum_b);
		comp_data_pos += leadNumArray_size;
		
		leadNumArray = comp_data_pos;
		convertByteArray2IntArray_fast_2b(reg_count, leadNumArray, leadNumArray_size, &leadNum_c);	
		comp_data_pos += leadNumArray_size;
				
		//reconstruct mid bytes...
		size_t * mid_byte_size_a = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* exactMidBytes_a = comp_data_pos;
		comp_data_pos += *mid_byte_size_a;
		
		size_t * mid_byte_size_b = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* exactMidBytes_b = comp_data_pos;	
		comp_data_pos += *mid_byte_size_b;	
		
		size_t * mid_byte_size_c = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* exactMidBytes_c = comp_data_pos;	
		comp_data_pos += *mid_byte_size_c;			
		
		//reconstruct the residualMidBits		
		size_t * resiMidBites_a_size = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* residualMidBits_a = comp_data_pos;
		comp_data_pos += *resiMidBites_a_size;

		size_t * resiMidBites_b_size = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* residualMidBits_b = comp_data_pos;
		comp_data_pos += *resiMidBites_b_size;
		
		size_t * resiMidBites_c_size = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* residualMidBits_c = comp_data_pos;
		comp_data_pos += *resiMidBites_c_size;
				
		//perform the decompression using the reconstructed leadNum, exactMidBytes and residualMidBits....
		decompressExactDataArray_float(leadNum_a, exactMidBytes_a, residualMidBits_a, reg_count, reqLength_a, medianValue_a, &dec_a);
		decompressExactDataArray_float(leadNum_b, exactMidBytes_b, residualMidBits_b, reg_count, reqLength_b, medianValue_b, &dec_b);
		decompressExactDataArray_float(leadNum_c, exactMidBytes_c, residualMidBits_c, reg_count, reqLength_c, medianValue_c, &dec_c);
		free(leadNum_a);
		free(leadNum_b);
		free(leadNum_c);
		
	}
	size_t total_unpred;
	memcpy(&total_unpred, comp_data_pos, sizeof(size_t));
	comp_data_pos += sizeof(size_t);
	float * unpred_data = (float *) comp_data_pos;
	comp_data_pos += total_unpred * sizeof(float);

	int * result_type = (int *) malloc(num_elements * sizeof(int));
	decode(comp_data_pos, num_elements, root, result_type);
	SZ_ReleaseHuffman(huffmanTree);
	
	int intvRadius = exe_params->intvRadius;
	
	int * type;

	float * data_pos = *data;
	size_t offset_x, offset_y;
	size_t current_blockcount_x, current_blockcount_y;
	size_t cur_unpred_count;

	float * dec_a_pos = dec_a;
	float * dec_b_pos = dec_b;
	float * dec_c_pos = dec_c;
	unsigned char * indicator_pos = indicator;
	if(use_mean){
		type = result_type;
		for(size_t i=0; i<num_x; i++){
			for(size_t j=0; j<num_y; j++){
				offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
				offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
				data_pos = *data + offset_x * dim0_offset + offset_y;

				current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
				current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;

				size_t current_block_elements = current_blockcount_x * current_blockcount_y;
				if(*indicator_pos){
					// decompress by SZ

					float * block_data_pos = data_pos;
					float pred;
					size_t index = 0;
					int type_;
					// d11 is current data
					size_t unpredictable_count = 0;
					float d00, d01, d10;
					for(size_t ii=0; ii<current_blockcount_x; ii++){
						for(size_t jj=0; jj<current_blockcount_y; jj++){
							type_ = type[index];
							if(type_ == intvRadius){
								*block_data_pos = mean;
							}
							else if(type_ == 0){
								*block_data_pos = unpred_data[unpredictable_count ++];
							}
							else{
								d00 = d01 = d10 = 1;
								if(i == 0 && ii == 0){
									d00 = d01 = 0;
								}
								if(j == 0 && jj == 0){
									d00 = d10 = 0;
								}
								if(d00){
									d00 = block_data_pos[- dim0_offset - 1];
								}
								if(d01){
									d01 = block_data_pos[- dim0_offset];
								}
								if(d10){
									d10 = block_data_pos[- 1];
								}
								if(type_ < intvRadius) type_ += 1;
								pred = d10 + d01 - d00;
								*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
							}
							index ++;
							block_data_pos ++;
						}
						block_data_pos += dim0_offset - current_blockcount_y;
					}
					cur_unpred_count = unpredictable_count;
				}
				else{
					// decompress by regression
					{
						float * block_data_pos = data_pos;
						float pred;
						int type_;
						size_t index = 0;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								type_ = type[index];
								if (type_ != 0){
									pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0];
									*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
								}
								else{
									*block_data_pos = unpred_data[unpredictable_count ++];
								}

								index ++;	
								block_data_pos ++;
							}
							block_data_pos += dim0_offset - current_blockcount_y;
						}
						cur_unpred_count = unpredictable_count;
					}

					dec_a_pos ++, dec_b_pos ++, dec_c_pos ++;
				}

				type += current_block_elements;
				indicator_pos ++;
				unpred_data += cur_unpred_count;
			}
		}
	}
	else{
		type = result_type;
		for(size_t i=0; i<num_x; i++){
			for(size_t j=0; j<num_y; j++){
				offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
				offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
				data_pos = *data + offset_x * dim0_offset + offset_y;

				current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
				current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;

				size_t current_block_elements = current_blockcount_x * current_blockcount_y;
				if(*indicator_pos){
					// decompress by SZ
					
					float * block_data_pos = data_pos;
					float pred;
					size_t index = 0;
					int type_;
					// d11 is current data
					size_t unpredictable_count = 0;
					float d00, d01, d10;
					for(size_t ii=0; ii<current_blockcount_x; ii++){
						for(size_t jj=0; jj<current_blockcount_y; jj++){
							type_ = type[index];
							if(type_ == 0){
								*block_data_pos = unpred_data[unpredictable_count ++];
							}
							else{
								d00 = d01 = d10 = 1;
								if(i == 0 && ii == 0){
									d00 = d01 = 0;
								}
								if(j == 0 && jj == 0){
									d00 = d10 = 0;
								}
								if(d00){
									d00 = block_data_pos[- dim0_offset - 1];
								}
								if(d01){
									d01 = block_data_pos[- dim0_offset];
								}
								if(d10){
									d10 = block_data_pos[- 1];
								}
								pred = d10 + d01 - d00;
								*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
							}
							index ++;
							block_data_pos ++;
						}
						block_data_pos += dim0_offset - current_blockcount_y;
					}
					cur_unpred_count = unpredictable_count;
				}
				else{
					// decompress by regression
					{
						float * block_data_pos = data_pos;
						float pred;
						int type_;
						size_t index = 0;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								type_ = type[index];
								if (type_ != 0){
									pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0];
									*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
								}
								else{
									*block_data_pos = unpred_data[unpredictable_count ++];
								}
								index ++;	
								block_data_pos ++;
							}
							block_data_pos += dim0_offset - current_blockcount_y;
						}
						cur_unpred_count = unpredictable_count;
					}
					dec_a_pos ++, dec_b_pos ++, dec_c_pos ++;
				}

				type += current_block_elements;
				indicator_pos ++;
				unpred_data += cur_unpred_count;
			}
		}
	}
	if(NULL != dec_a) free(dec_a);
	if(NULL != dec_b) free(dec_b);
	if(NULL != dec_c) free(dec_c);

	free(indicator);
	free(result_type);
}


void decompressDataSeries_float_3D_nonblocked_with_blocked_regression(float** data, size_t r1, size_t r2, size_t r3, unsigned char* comp_data){

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	size_t num_elements = r1 * r2 * r3;

	*data = (float*)malloc(sizeof(float)*num_elements);

	unsigned char * comp_data_pos = comp_data;

	size_t block_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += sizeof(int);
	// calculate block dims
	size_t num_x, num_y, num_z;
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r2, num_y, block_size);
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r3, num_z, block_size);

	size_t split_index_x, split_index_y, split_index_z;
	size_t early_blockcount_x, early_blockcount_y, early_blockcount_z;
	size_t late_blockcount_x, late_blockcount_y, late_blockcount_z;
	SZ_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	SZ_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);
	SZ_COMPUTE_BLOCKCOUNT(r3, num_z, split_index_z, early_blockcount_z, late_blockcount_z);

	size_t num_blocks = num_x * num_y * num_z;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += sizeof(double);
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += sizeof(int);

	updateQuantizationInfo(intervals);

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += sizeof(int);
	
	int stateNum = 2*intervals;
	HuffmanTree* huffmanTree = createHuffmanTree(stateNum);	
	
	int nodeCount = bytesToInt_bigEndian(comp_data_pos);
	node root = reconstruct_HuffTree_from_bytes_anyStates(huffmanTree,comp_data_pos+4, nodeCount);
	comp_data_pos += sizeof(int) + tree_size;

	float mean;
	unsigned char use_mean;
	memcpy(&use_mean, comp_data_pos, sizeof(unsigned char));
	comp_data_pos += sizeof(unsigned char);
	memcpy(&mean, comp_data_pos, sizeof(float));
	comp_data_pos += sizeof(float);
	size_t reg_count = 0;

	unsigned char * indicator;
	size_t indicator_bitlength = (num_blocks - 1)/8 + 1;
	convertByteArray2IntArray_fast_1b(num_blocks, comp_data_pos, indicator_bitlength, &indicator);
	comp_data_pos += indicator_bitlength;
	for(size_t i=0; i<num_blocks; i++){
		if(!indicator[i]) reg_count ++;
	}
	float* dec_a = NULL, *dec_b = NULL, *dec_c = NULL, *dec_d = NULL;
	if(reg_count > 0){
		float * medians = (float *) comp_data_pos;
		float medianValue_a = *medians;
		float medianValue_b = *(medians + 1);
		float medianValue_c = *(medians + 2);
		float medianValue_d = *(medians + 3);
		comp_data_pos += 4 * sizeof(float);

		int * reqLength = (int *) comp_data_pos;
		int reqLength_a = *reqLength;
		int reqLength_b = *(reqLength + 1);
		int reqLength_c = *(reqLength + 2);
		int reqLength_d = *(reqLength + 3);
		comp_data_pos += 4 * sizeof(int);
		
		//reconstruct leading_number array in the form of meaningful integers....
		size_t leadNumArray_size = (reg_count-1)/4+1;
		
		unsigned char* leadNum_a = NULL, *leadNum_b = NULL, *leadNum_c = NULL, *leadNum_d = NULL;
		
		unsigned char* leadNumArray = comp_data_pos;
		convertByteArray2IntArray_fast_2b(reg_count, leadNumArray, leadNumArray_size, &leadNum_a);
		comp_data_pos += leadNumArray_size;
		
		leadNumArray = comp_data_pos;
		convertByteArray2IntArray_fast_2b(reg_count, leadNumArray, leadNumArray_size, &leadNum_b);
		comp_data_pos += leadNumArray_size;
		
		leadNumArray = comp_data_pos;
		convertByteArray2IntArray_fast_2b(reg_count, leadNumArray, leadNumArray_size, &leadNum_c);	
		comp_data_pos += leadNumArray_size;
		
		leadNumArray = comp_data_pos;
		convertByteArray2IntArray_fast_2b(reg_count, leadNumArray, leadNumArray_size, &leadNum_d);			
		comp_data_pos += leadNumArray_size;
		
		//reconstruct mid bytes...
		size_t * mid_byte_size_a = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* exactMidBytes_a = comp_data_pos;
		comp_data_pos += *mid_byte_size_a;
		
		size_t * mid_byte_size_b = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* exactMidBytes_b = comp_data_pos;	
		comp_data_pos += *mid_byte_size_b;	
		
		size_t * mid_byte_size_c = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* exactMidBytes_c = comp_data_pos;	
		comp_data_pos += *mid_byte_size_c;			
		
		size_t * mid_byte_size_d = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* exactMidBytes_d = comp_data_pos;	
		comp_data_pos += *mid_byte_size_d;			

		//reconstruct the residualMidBits		
		size_t * resiMidBites_a_size = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* residualMidBits_a = comp_data_pos;
		comp_data_pos += *resiMidBites_a_size;

		size_t * resiMidBites_b_size = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* residualMidBits_b = comp_data_pos;
		comp_data_pos += *resiMidBites_b_size;
		
		size_t * resiMidBites_c_size = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* residualMidBits_c = comp_data_pos;
		comp_data_pos += *resiMidBites_c_size;
		
		size_t * resiMidBites_d_size = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* residualMidBits_d = comp_data_pos;
		comp_data_pos += *resiMidBites_d_size;				
		
		//perform the decompression using the reconstructed leadNum, exactMidBytes and residualMidBits....
		decompressExactDataArray_float(leadNum_a, exactMidBytes_a, residualMidBits_a, reg_count, reqLength_a, medianValue_a, &dec_a);
		decompressExactDataArray_float(leadNum_b, exactMidBytes_b, residualMidBits_b, reg_count, reqLength_b, medianValue_b, &dec_b);
		decompressExactDataArray_float(leadNum_c, exactMidBytes_c, residualMidBits_c, reg_count, reqLength_c, medianValue_c, &dec_c);
		decompressExactDataArray_float(leadNum_d, exactMidBytes_d, residualMidBits_d, reg_count, reqLength_d, medianValue_d, &dec_d);

		free(leadNum_a);
		free(leadNum_b);
		free(leadNum_c);
		free(leadNum_d);				
		
	}

	size_t total_unpred;
	memcpy(&total_unpred, comp_data_pos, sizeof(size_t));
	comp_data_pos += sizeof(size_t);
	float * unpred_data = (float *) comp_data_pos;
	comp_data_pos += total_unpred * sizeof(float);

	int * result_type = (int *) malloc(num_elements * sizeof(int));
	decode(comp_data_pos, num_elements, root, result_type);
	SZ_ReleaseHuffman(huffmanTree);
	
	int intvRadius = exe_params->intvRadius;
	
	int * type;
	float * data_pos = *data;
	size_t offset_x, offset_y, offset_z;
	size_t current_blockcount_x, current_blockcount_y, current_blockcount_z;
	size_t cur_unpred_count;
	float * dec_a_pos = dec_a;
	float * dec_b_pos = dec_b;
	float * dec_c_pos = dec_c;
	float * dec_d_pos = dec_d;
	unsigned char * indicator_pos = indicator;
	if(use_mean){
		type = result_type;
		// i == 0
		{
			// j == 0
			{
				// k == 0
				{
					data_pos = *data;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							// jj == 0
							{
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = 0;
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;						
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim0_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				// i == 0 j == 0 k != 0
				for(size_t k=1; k<num_z; k++){
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_z;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							// jj == 0
							{
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}// end j==0
			for(size_t j=1; j<num_y; j++){
				// k == 0
				{
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					data_pos = *data + offset_y * dim1_offset;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				for(size_t k=1; k<num_z; k++){
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_y * dim1_offset + offset_z;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}
		} // end i==0
		for(size_t i=1; i<num_x; i++){
			// j == 0
			{
				// k == 0
				{
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					data_pos = *data + offset_x * dim0_offset;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim0_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				for(size_t k=1; k<num_z; k++){
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_x * dim0_offset + offset_z;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}// end j = 0
			for(size_t j=1; j<num_y; j++){
				// k == 0
				{
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					data_pos = *data + offset_x * dim0_offset + offset_y * dim1_offset;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						// reg_params += 4;
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				for(size_t k=1; k<num_z; k++){
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}
		}
	}
	else{
		type = result_type;
		// i == 0
		{
			// j == 0
			{
				// k == 0
				{
					data_pos = *data;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							// jj == 0
							{
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = 0;
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;						
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim0_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				// i == 0 j == 0 k != 0
				for(size_t k=1; k<num_z; k++){
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_z;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							// jj == 0
							{
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}// end j==0
			for(size_t j=1; j<num_y; j++){
				// k == 0
				{
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					data_pos = *data + offset_y * dim1_offset;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				for(size_t k=1; k<num_z; k++){
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_y * dim1_offset + offset_z;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}
		} // end i==0
		for(size_t i=1; i<num_x; i++){
			// j == 0
			{
				// k == 0
				{
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					data_pos = *data + offset_x * dim0_offset;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim0_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				for(size_t k=1; k<num_z; k++){
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_x * dim0_offset + offset_z;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}// end j = 0
			for(size_t j=1; j<num_y; j++){
				// k == 0
				{
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					data_pos = *data + offset_x * dim0_offset + offset_y * dim1_offset;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				for(size_t k=1; k<num_z; k++){
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}
		}
	}
	if(NULL != dec_a) free(dec_a);
	if(NULL != dec_b) free(dec_b);
	if(NULL != dec_c) free(dec_c);
	if(NULL != dec_d) free(dec_d);

	free(indicator);
	free(result_type);
}
