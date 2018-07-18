/**
 *  @file test_compress_ts.c
 *  @author Sheng Di
 *  @date May, 2018
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include "sz.h"
#include "rw.h"

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;

#define NB_variable 6


void cost_start()
{
	totalCost = 0;
        gettimeofday(&costStart, NULL);
}

void cost_end()
{
        double elapsed;
        struct timeval costEnd;
        gettimeofday(&costEnd, NULL);
        elapsed = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
        totalCost += elapsed;
}


int main(int argc, char * argv[])
{
    int i = 0;
    size_t r5=0,r4=0,r3=0,r2=0,r1=0;
    char oriDir[640], outputDir[640], outputFilePath[600];
    char *cfgFile;
    float eb = 0.001;
    
    if(argc < 3)
    {
		printf("Test case: testfloat_compress_ts [config_file] [srcDir] [dimension sizes...]\n");
		printf("Example: testfloat_compress_ts sz.config /home/sdi/Data/Hurricane-ISA/consecutive-steps 500 500 100\n");
		exit(0);
    }
   
    cfgFile=argv[1];
    sprintf(oriDir, "%s", argv[2]);//add a comment to test github
    if(argc>=4)
		r1 = atoi(argv[3]); //8
    if(argc>=5)
		//r2 = atoi(argv[4]); //8
        eb = atof(argv[4]);
    if(argc>=6)
		r3 = atoi(argv[5]); //128
    if(argc>=7)
        r4 = atoi(argv[6]);
    if(argc>=8)
        r5 = atoi(argv[7]);
   
    printf("cfgFile=%s\n", cfgFile); 
    int status = SZ_Init(cfgFile);
    if(status == SZ_NSCS)
		exit(0);
    sprintf(outputDir, "%s", oriDir);
   
    char oriFilePath[600];
    size_t nbEle;
    size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
    //float *data = (float*)malloc(sizeof(float)*dataLength);
    float **data = (float**) malloc(NB_variable * sizeof(float*));
    for (i = 0; i < NB_variable; i++){
        data[i] = (float*) malloc(dataLength * sizeof(float));
    }
    int64_t* index = (int64_t*) malloc(dataLength * sizeof(int64_t));
    //SZ_registerVar("CLOUDf", SZ_FLOAT, data, REL, 0, 0.001, 0, r5, r4, r3, r2, r1);

    SZ_registerVar("x", SZ_FLOAT, data[0], REL, 0, eb, 0, r5, r4, r3, r2, r1);
    SZ_registerVar("y", SZ_FLOAT, data[1], REL, 0, eb, 0, r5, r4, r3, r2, r1);
    SZ_registerVar("z", SZ_FLOAT, data[2], REL, 0, eb, 0, r5, r4, r3, r2, r1);
    SZ_registerVar("vx", SZ_FLOAT, data[3], REL, 0, eb, 0, r5, r4, r3, r2, r1);
    SZ_registerVar("vy", SZ_FLOAT, data[4], REL, 0, eb, 0, r5, r4, r3, r2, r1);
    SZ_registerVar("vz", SZ_FLOAT, data[5], REL, 0, eb, 0, r5, r4, r3, r2, r1);
    SZ_registerVar("index", SZ_INT64, index, REL, 0, eb, 0, r5, r4, r3, r2, r1);

    if(status != SZ_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }

    int file_num[6] = {100, 102, 105, 107, 110, 113};
   
    size_t outSize; 
    unsigned char *bytes = NULL;
    for(i=0;i<6;i++)
	{
		printf("simulation time step %d\n", i);
        int m = 0;
        for (m = 0; m < NB_variable; m++){
            sprintf(oriFilePath, "%s/m000.full.mpicosmo.%d#21-%d.dat", oriDir, file_num[i], m);
            float* data_ = readFloatData(oriFilePath, &nbEle, &status);
            memcpy(data[m], data_, nbEle*sizeof(float));
            free(data_);
        }
        sprintf(oriFilePath, "%s/m000.full.mpicosmo.%d#21-i.dat", oriDir, file_num[i]);
        int64_t* index_ = readInt64Data(oriFilePath, &nbEle, &status);
        memcpy(index, index_, nbEle*sizeof(int64_t));
        //printf("--------------- %lld \n", index_[0]);
		//float *data_ = readFloatData(oriFilePath, &nbEle, &status);
		//memcpy(data, data_, nbEle*sizeof(float));
		cost_start();
		SZ_compress_ts_vlct(&bytes, &outSize);
		cost_end();
		printf("timecost=%f\n",totalCost); 
		sprintf(outputFilePath, "%s/QCLOUDf%02d.bin.dat.sz2", outputDir, i);
		printf("writing compressed data to %s\n", outputFilePath);
		writeByteData(bytes, outSize, outputFilePath, &status); 
		free(bytes);
		//free(data_);
	}
    
    printf("done\n");
    free(data);
    SZ_Finalize();
    
    return 0;
}
