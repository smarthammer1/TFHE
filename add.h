#include <tfhe/numeric_functions.h>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_core.h>
#include <tfhe/tfhe_io.h>
#include <math.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/lweparams.h>
#include <tfhe/lwekey.h>
#include <tfhe/lwesamples.h>
#include <tfhe/lwekeyswitch.h>
#include <tfhe/lwe-functions.h>
#include <tfhe/lwebootstrappingkey.h>
#include <iostream>
#include <vector>
#include <thread>

using namespace std;
int32_t get_N();

// void HomAdd(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk);

void tfhe_bootstrap_rotated_woKS_FFT(LweSample *result,
                                    const LweBootstrappingKeyFFT *bk,
                                    Torus32 mu,
                                    const LweSample *x);

void tfhe_rotate_bootstrap_FFT(LweSample *result,
                               const LweBootstrappingKeyFFT *bk,
                               Torus32 mu,
                               const LweSample *x);

void bootsCarry(LweSample* result, const LweSample* ca, const LweSample* cb, const LweSample* cc, const TFheGateBootstrappingCloudKeySet *bk);
void bootsSum(LweSample* result, const LweSample* ca, const LweSample* cb, const LweSample* cc, const TFheGateBootstrappingCloudKeySet *bk);
void newADD(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk);

void bootsBorrow(LweSample* result, const LweSample* ca, const LweSample* cb, const LweSample* cc, const TFheGateBootstrappingCloudKeySet *bk);
void bootsCompS(LweSample* result, const LweSample* ca, const LweSample* cb, const LweSample* cc, const TFheGateBootstrappingCloudKeySet *bk);
void bootsCompB(LweSample* result, const LweSample* ca, const LweSample* cb, const LweSample* cc, const TFheGateBootstrappingCloudKeySet *bk);
// void BootsSort(LweSample* res1, LweSample* res2, LweSample* a, LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk);

void newSUB(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newABS(LweSample* res, const LweSample* a, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newCompS(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newCompB(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newCompSE(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newCompLE(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newMax(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newMin(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk);


void newMultiReal(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newMultiReal2(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newRealDiv(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newMult(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk);

void HomCovar(LweSample* res, LweSample** a, LweSample** b, const int length, const int N, const TFheGateBootstrappingCloudKeySet* bk);
void newCovar(LweSample* res, LweSample** a, LweSample** b, const int length, const int N, const TFheGateBootstrappingCloudKeySet* bk);
double cov(double *plain_a, double *plain_b, int32_t N);

void newMaxValue(LweSample* res, vector<vector<LweSample*>> matrix, const int length, const TFheGateBootstrappingCloudKeySet* bk);