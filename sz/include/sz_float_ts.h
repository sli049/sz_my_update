/**
 *  @file sz_float_ts.h
 *  @author Sheng Di
 *  @date May, 2018
 *  @brief Header file for the sz_float_ts.c.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
#include "TightDataPointStorageF.h"

#ifndef _SZ_Float_TS_H
#define _SZ_Float_TS_H

#ifdef __cplusplus
extern "C" {
#endif
unsigned int optimize_intervals_float_1D_ts(float *oriData, size_t dataLength, float* preData, double realPrecision);

TightDataPointStorageF* SZ_compress_float_1D_MDQ_ts(float *oriData, size_t dataLength, sz_multisteps* multisteps,
double realPrecision, float valueRangeSize, float medianValue_f);

TightDataPointStorageF* SZ_compress_float_1D_MDQ_ts_vlct(float *oriData, size_t dataLength, sz_multisteps* multisteps,
double realPrecision, float valueRangeSize, float medianValue_f);//sihuan added, in this, v is used to predict x more accurately

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _SZ_Float_TS_H  ----- */

