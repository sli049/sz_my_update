/**
 *  @file szd_float_ts.c
 *  @author Sheng Di and Dingwen Tao
 *  @date Aug, 2016
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
#include "szd_float_ts.h"

void decompressDataSeries_float_1D_ts(float** data, size_t dataSeriesLength, sz_multisteps* multisteps, TightDataPointStorageF* tdps) 
{
	float* lastSnapshotData = (float*)multisteps->hist_data;
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
	float medianValue, exactData, predValue = 0;
	
	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;

	//now decompress the bitarray and use it to do decompression(this should be optimized to execute for one step while not one varialble);
	char bit_file_name[50];
	sprintf(bit_file_name, "bitarray_out_sz_%d.sz1", sz_tsc->currentStep);
	size_t byteLength_tmp;
	int status_tmp;
	unsigned char* bit_read = readByteData(bit_file_name, &byteLength_tmp, &status_tmp);
	size_t oriSize = bytesToSize(bit_read);
	bit_read += sizeof(size_t);
	unsigned long cmpSize = (unsigned long) bytesToLong_bigEndian(bit_read);;
	bit_read += sizeof(unsigned long);

	if(cmpSize == byteLength_tmp - sizeof(size_t) - sizeof(unsigned long))
		printf("the two zlib compuression size match\n");
	else("the two zlib compression size DON'T match\n");

	unsigned char* dec_bit_arr = (unsigned char*) malloc(sizeof(unsigned char) * oriSize);
	zlib_uncompress(bit_read, cmpSize, &dec_bit_arr, oriSize);


	int m = 0;
	while (dec_bit_arr[m] == '1'){m++;}
	
	int type_;
	for (i = 0; i < dataSeriesLength; i++, m++) { //need to j++ also
		type_ = type[i];
		while(dec_bit_arr[m] == '1'){m++;}
		switch (type_) {
		case 0:

			//increase the j pointer
			
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
			//predValue = (*data)[i-1];
			if(confparams_dec->szMode == SZ_TEMPORAL_COMPRESSION)
				//predValue = lastSnapshotData[i];//change the i pointer to be j here
				predValue = lastSnapshotData[m];
			(*data)[i] = predValue + (type_-exe_params->intvRadius)*interval;
			break;
		}
		//printf("%.30G\n",(*data)[i]);
	}
	
	memcpy(multisteps->hist_data, (*data), dataSeriesLength*sizeof(float));
	
	free(leadNum);
	free(type);
	return;
}


void decompressDataSeries_float_1D_ts_vlct(float** data, size_t dataSeriesLength, sz_multisteps* multisteps, TightDataPointStorageF* tdps) 
{
	float* lastSnapshotData = (float*)multisteps->hist_data;

	SZ_Variable* v_tmp = NULL;
	v_tmp = sz_varset->header->next;
	if (!strcmp(v_global->varName, "x")){
		printf("vlct find vx for x\n");
		while(strcmp(v_tmp->varName, "vx")) {v_tmp = v_tmp->next;}
		printf("found varible name: %s\n", v_tmp->varName);
	}
	else if (!strcmp(v_global->varName, "y")){
		printf("vlct find vy for y\n");
		while(strcmp(v_tmp->varName, "vy")) {v_tmp = v_tmp->next;}
		printf("found varible name: %s\n", v_tmp->varName);
	}
	else if (!strcmp(v_global->varName, "z")){
		printf("vlct find vz for z\n");
		while(strcmp(v_tmp->varName, "vz")) {v_tmp = v_tmp->next;}
		printf("found varible name: %s\n", v_tmp->varName);
	}
	else printf("!!!!!!!!Should not call SZ_compress_float_1D_MDQ_ts_vlct, instead consider SZ_compress_float_1D_MDQ_ts!!!!!\n");

	float* lastSnapshotData_vol = (float*)(v_tmp->multisteps->hist_data);



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
	float medianValue, exactData, predValue = 0;
	
	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;

	//now decompress the bitarray and use it to do decompression(this should be optimized to execute for one step while not one varialble);
	char bit_file_name[50];
	sprintf(bit_file_name, "bitarray_out_sz_%d.sz1", sz_tsc->currentStep);
	size_t byteLength_tmp;
	int status_tmp;
	unsigned char* bit_read = readByteData(bit_file_name, &byteLength_tmp, &status_tmp);
	size_t oriSize = bytesToSize(bit_read);
	bit_read += sizeof(size_t);
	unsigned long cmpSize = (unsigned long) bytesToLong_bigEndian(bit_read);;
	bit_read += sizeof(unsigned long);

	if(cmpSize == byteLength_tmp - sizeof(size_t) - sizeof(unsigned long))
		printf("the two zlib compuression size match\n");
	else("the two zlib compression size DON'T match\n");

	unsigned char* dec_bit_arr = (unsigned char*) malloc(sizeof(unsigned char) * oriSize);
	zlib_uncompress(bit_read, cmpSize, &dec_bit_arr, oriSize);


	int m = 0;
	while (dec_bit_arr[m] == '1'){m++;}
	
	int type_;
	float delta_t = delta_t_opt[sz_tsc->currentStep-1];
	for (i = 0; i < dataSeriesLength; i++, m++) { //need to j++ also
		type_ = type[i];
		while(dec_bit_arr[m] == '1'){m++;}
		switch (type_) {
		case 0:

			//increase the j pointer
			
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
			//predValue = (*data)[i-1];
			if(confparams_dec->szMode == SZ_TEMPORAL_COMPRESSION)
				//predValue = lastSnapshotData[i];//change the i pointer to be j here
				predValue = lastSnapshotData[m] + lastSnapshotData_vol[m]*delta_t;
			(*data)[i] = predValue + (type_-exe_params->intvRadius)*interval;
			break;
		}
		//printf("%.30G\n",(*data)[i]);
	}
	
	memcpy(multisteps->hist_data, (*data), dataSeriesLength*sizeof(float));
	
	free(leadNum);
	free(type);
	return;
}
