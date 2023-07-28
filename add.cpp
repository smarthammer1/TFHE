#include "add.h"
#include "Comparator.h"
#include <tfhe/tfhe_gate_bootstrapping_functions.h>
#include "HomOper.c"
#include "newfile.h"

int32_t get_N();

void tfhe_bootstrap_rotated_woKS_FFT(LweSample *result,
                                    const LweBootstrappingKeyFFT *bk,
                                    Torus32 mu,
                                    const LweSample *x) {

    const TGswParams *bk_params = bk->bk_params;
    const TLweParams *accum_params = bk->accum_params;
    const LweParams *in_params = bk->in_out_params;
    const int32_t N = accum_params->N;
    const int32_t Nx2 = 2 * N;
    const int32_t n = in_params->n;
    const int32_t halfN = N/2;

    TorusPolynomial *testvect = new_TorusPolynomial(N);
    int32_t *bara = new int32_t[N];


    // Modulus switching
    int32_t barb = modSwitchFromTorus32(x->b, Nx2);
    for (int32_t i = 0; i < n; i++) {
        bara[i] = modSwitchFromTorus32(x->a[i], Nx2);
    }

    for (int32_t i = 0; i < halfN; i++) testvect->coefsT[i] = -mu;
    for (int32_t i = halfN; i < N; i++) testvect->coefsT[i] = mu;
    

    // Bootstrapping rotation and extraction
    tfhe_blindRotateAndExtract_FFT(result, testvect, bk->bkFFT, barb, bara, n, bk_params);


    delete[] bara;
    delete_TorusPolynomial(testvect);
}

void tfhe_rotate_bootstrap_FFT(LweSample *result,
                               const LweBootstrappingKeyFFT *bk,
                               Torus32 mu,
                               const LweSample *x) {

    LweSample *u = new_LweSample(&bk->accum_params->extracted_lweparams);

    tfhe_bootstrap_rotated_woKS_FFT(u, bk, mu, x);
    // Key switching
    lweKeySwitch(result, bk->ks, u);

    delete_LweSample(u);
}

void bootsCarry(LweSample* result, const LweSample* ca, const LweSample* cb, const LweSample* cc, const TFheGateBootstrappingCloudKeySet *bk) {
    static const Torus32 MU = modSwitchToTorus32(1, 8);
    const LweParams *in_out_params = bk->params->in_out_params;

    LweSample *temp_result = new_LweSample(in_out_params);

    // static const Torus32 addtempConst = modSwitchToTorus32(2, 8);
    // lweNoiselessTrivial(temp_result, addtempConst, in_out_params);
    lweClear(temp_result, in_out_params);
    lweAddTo(temp_result, ca, in_out_params);
    lweAddTo(temp_result, cb, in_out_params);
    lweAddTo(temp_result, cc, in_out_params);

    // tfhe_bootstrap_rotated_woKS_FFT(result, bk->bkFFT, MU, temp_result);
    tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);

    
    delete_LweSample(temp_result);
}

void bootsSum(LweSample* result, const LweSample* ca, const LweSample* cb, const LweSample* cc, const TFheGateBootstrappingCloudKeySet *bk) {
    static const Torus32 MU = modSwitchToTorus32(1, 8);
    const LweParams *in_out_params = bk->params->in_out_params;

    LweSample *temp_result = new_LweSample(in_out_params);

    static const Torus32 addresConst = modSwitchToTorus32(1, 2);
    lweNoiselessTrivial(temp_result, addresConst, in_out_params);
    lweAddMulTo(temp_result, 2, ca, in_out_params);
    lweAddMulTo(temp_result, 2, cb, in_out_params);
    lweAddMulTo(temp_result, 2, cc, in_out_params);

    tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);
    // tfhe_rotate_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);
    
    delete_LweSample(temp_result);
}

void newADD(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
	LweSample* temp = new_LweSample(bk->params->in_out_params);
    LweSample* tmp_a = new_LweSample(bk->params->in_out_params); 
    LweSample* tmp_b = new_LweSample(bk->params->in_out_params); 
    bootsCONSTANT(temp, 0, bk);   

    // bootsXOR(&res[0], &a[0], &b[0], bk);
    // bootsAND(temp, &a[0], &b[0], bk);
    for(int i = 0; i < length-1; i++){
        bootsCOPY(tmp_a, &a[i], bk);
        bootsCOPY(tmp_b, &b[i], bk);
        bootsSum(&res[i], tmp_a, tmp_b, temp, bk);
        bootsCarry(temp, tmp_a, tmp_b, temp, bk);       //fixed _  putting res in input a or b again : all ok!
    }
    bootsSum(&res[length-1], &a[length-1], &b[length-1], temp, bk);
    delete_LweSample(temp);
    delete_LweSample(tmp_a);
    delete_LweSample(tmp_b);
}

void bootsBorrow(LweSample* result, const LweSample* ca, const LweSample* cb, const LweSample* cc, const TFheGateBootstrappingCloudKeySet *bk) {
    static const Torus32 MU = modSwitchToTorus32(1, 8);
    const LweParams *in_out_params = bk->params->in_out_params;

    LweSample *temp_result = new_LweSample(in_out_params);

    lweClear(temp_result, in_out_params);
    lweSubTo(temp_result, ca, in_out_params);
    lweAddTo(temp_result, cb, in_out_params);
    lweAddTo(temp_result, cc, in_out_params);

    tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);
    delete_LweSample(temp_result);
}

// EXACTLY same as bootsBorrow 
void bootsCompS(LweSample* result, const LweSample* ca, const LweSample* cb, const LweSample* cc, const TFheGateBootstrappingCloudKeySet *bk) {
    static const Torus32 MU = modSwitchToTorus32(1, 8);
    const LweParams *in_out_params = bk->params->in_out_params;

    LweSample *temp_result = new_LweSample(in_out_params);

    lweClear(temp_result, in_out_params);
    lweSubTo(temp_result, ca, in_out_params);
    lweAddTo(temp_result, cb, in_out_params);
    lweAddTo(temp_result, cc, in_out_params);

    tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);
    delete_LweSample(temp_result);
}

void bootsCompB(LweSample* result, const LweSample* ca, const LweSample* cb, const LweSample* cc, const TFheGateBootstrappingCloudKeySet *bk) {
    static const Torus32 MU = modSwitchToTorus32(1, 8);
    const LweParams *in_out_params = bk->params->in_out_params;

    LweSample *temp_result = new_LweSample(in_out_params);

    lweClear(temp_result, in_out_params);
    lweAddTo(temp_result, ca, in_out_params);
    lweSubTo(temp_result, cb, in_out_params);
    lweAddTo(temp_result, cc, in_out_params);

    tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);
    delete_LweSample(temp_result);
}

void newSUB(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
	LweSample* temp = new_LweSample(bk->params->in_out_params);	 
    LweSample* tmp = new_LweSample(bk->params->in_out_params);
    bootsCONSTANT(temp, 0, bk);   

    // bootsXOR(&res[0], &a[0], &b[0], bk);
    // bootsAND(temp, &a[0], &b[0], bk);

    for(int i = 0; i < length-1; i++){
        bootsCOPY(tmp, &a[i], bk);
        bootsSum(&res[i], tmp, &b[i], temp, bk);
        bootsBorrow(temp, tmp, &b[i], temp, bk);
    }
    bootsSum(&res[length-1], &a[length-1], &b[length-1], temp, bk);
    delete_LweSample(temp);
    delete_LweSample(tmp);
}

void newCompS(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    LweSample* tmp = new_gate_bootstrapping_ciphertext_array(length, bk->params);

    newSUB(tmp, a, b, length, bk);
    bootsCOPY(res, &tmp[length-1], bk);
    delete_gate_bootstrapping_ciphertext(tmp);
}

void newCompB(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    LweSample* tmp = new_gate_bootstrapping_ciphertext_array(length, bk->params);

    newSUB(tmp, b, a, length, bk);
    bootsCOPY(res, &tmp[length-1], bk);
    delete_gate_bootstrapping_ciphertext(tmp);
}

void newCompSE(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    
    LweSample* tmp = new_gate_bootstrapping_ciphertext(bk->params);
    newCompB(tmp, a, b, length, bk);
    bootsNOT(res, tmp, bk);

    delete_gate_bootstrapping_ciphertext(tmp);
}

void newCompLE(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    
    LweSample* tmp = new_gate_bootstrapping_ciphertext(bk->params);
    newCompS(tmp, a, b, length, bk);
    bootsNOT(res, tmp, bk);

    delete_gate_bootstrapping_ciphertext(tmp);
}

void newMax(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk){
    
    LweSample* tmp = new_gate_bootstrapping_ciphertext(bk->params);

    newCompB(tmp, a, b, length, bk);

    for (int i = 0; i < length; i++) {
        bootsMUX(&res[i], tmp, &a[i], &b[i], bk);
    }
    delete_gate_bootstrapping_ciphertext(tmp);
}

void newReLUMax(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    
    LweSample* tmp = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* sign = new_gate_bootstrapping_ciphertext(bk->params);

    newSUB(tmp, a, b, length, bk);
    bootsNOT(sign, &tmp[length-1], bk);
    for(int i = 0; i < length; i++) {
        bootsAND(&tmp[i], &tmp[i], sign, bk);
    }
    newADD(res, tmp, b, length, bk);

    delete_gate_bootstrapping_ciphertext(sign);
    delete_gate_bootstrapping_ciphertext_array(length, tmp);

}

void newMin(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk){
    
    LweSample* tmp = new_gate_bootstrapping_ciphertext(bk->params);

    newCompS(tmp, a, b, length, bk);

    for (int i = 0; i < length; i++) {
        bootsMUX(&res[i], tmp, &a[i], &b[i], bk);
    }
    delete_gate_bootstrapping_ciphertext(tmp);
}

void newABS(LweSample* res, const LweSample* a, const int length, const TFheGateBootstrappingCloudKeySet* bk) {

    LweSample* tmp = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    for(int i = 1; i < length; i++) {
        bootsCONSTANT(&tmp[i], 0, bk);
    }

    bootsCOPY(&tmp[0], &a[length-1], bk);

    for(int i = 0; i < length; i++) {
        bootsXOR(&res[i], &a[i], &tmp[0], bk);
    }

    newADD(res, res, tmp, length, bk);

    delete_gate_bootstrapping_ciphertext(tmp);
}

void newMultiReal(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
	
	LweSample* temp = new_gate_bootstrapping_ciphertext_array(2, bk->params);
	LweSample* A = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* AA = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* B = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* C = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* D = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* E = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* F = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	
	newABS(A, a, length, bk);
	newABS(B, b, length, bk);
	for(int i = 0; i < length; i++){
		bootsAND(&AA[i], &A[i], &B[0], bk);}
	HomLShift(AA, AA, length, length/2, bk);	


	
	for(int i = 1; i < length-1; i++){
		if(i < length/2){
			HomLShift(C, A, length, length/2-i, bk);

			for(int j = 0; j < length; j++){
				bootsAND(&D[j], &C[j], &B[i], bk);
            }
			newADD(AA, AA, D, length, bk);
		}
		else if(i == length/2){
			for(int j = 0; j < length; j++)
				bootsAND(&D[j], &A[j], &B[i], bk);

			newADD(AA, AA, D, length, bk);
		}

		else {
			HomRShift(C, A, length, i-length/2, bk);
			for(int j = 0; j < length; j++){
				bootsAND(&D[j], &C[j], &B[i], bk);}

			newADD(AA, AA, D, length, bk);
		}

	}
	bootsCONSTANT(&AA[length-1], 0, bk);

	HomTwosComplement(D, AA, length, bk);
	
	bootsXOR(&temp[0], &a[length-1], &b[length-1], bk);
	bootsNOT(&temp[1], &temp[0], bk);

	for(int i = 0; i < length; i++){
		bootsAND(&E[i], &AA[i], &temp[1], bk);
		bootsAND(&F[i], &D[i], &temp[0], bk);}

	newADD(res, E, F, length, bk);

	delete_gate_bootstrapping_ciphertext_array(2, temp);
	delete_gate_bootstrapping_ciphertext_array(length, A);
	delete_gate_bootstrapping_ciphertext_array(length, AA);
	delete_gate_bootstrapping_ciphertext_array(length, B);
	delete_gate_bootstrapping_ciphertext_array(length, C);
	delete_gate_bootstrapping_ciphertext_array(length, D);
	delete_gate_bootstrapping_ciphertext_array(length, E);
	delete_gate_bootstrapping_ciphertext_array(length, F);
}


void newMultiReal2(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
	
	LweSample* temp = new_gate_bootstrapping_ciphertext_array(2, bk->params);
	LweSample* A = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* AA = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* B = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* C = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* D = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* E = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* F = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	
	newABS(A, a, length, bk);
	newABS(B, b, length, bk);
    for(int i = 0; i < length; i++){
        bootsCONSTANT(&AA[i], 0, bk);
        bootsCONSTANT(&D[i], 0, bk);
    }

	for(int i = 0; i < length/2 - 1; i++){
		bootsAND(&AA[i], &A[i + length/2], &B[0], bk);
    }
	// HomLShift(AA, AA, length, length/2, bk);	


	
	for(int i = 1; i < length-1; i++){
		if(i < length/2){

			for(int j = 0; j < length/2 - 1 + i; j++){
				bootsAND(&D[j], &A[length/2 - i + j], &B[i], bk);
            }
			newADD(AA, AA, D, length/2 + i, bk);
		}
		else if(i == length/2){
			for(int j = 0; j < length-1; j++){
				bootsAND(&D[j], &A[j], &B[i], bk);
            }
			newADD(AA, AA, D, length, bk);
		}

		else {
			for(int j = 0; j < 3*(length/2) - i - 1; j++){
				bootsAND(&D[j], &A[j], &B[i], bk);
            }

			newADD(AA + (i - length/2), AA + (i - length/2), D, 3*(length/2) - i - 1, bk);
		}

	}
	bootsCONSTANT(&AA[length-1], 0, bk);

	HomTwosComplement(D, AA, length, bk);
	
	bootsXOR(&temp[0], &a[length-1], &b[length-1], bk);
	bootsNOT(&temp[1], &temp[0], bk);

	for(int i = 0; i < length; i++){
		bootsAND(&E[i], &AA[i], &temp[1], bk);
		bootsAND(&F[i], &D[i], &temp[0], bk);}

	newADD(res, E, F, length, bk);

	delete_gate_bootstrapping_ciphertext_array(2, temp);
	delete_gate_bootstrapping_ciphertext_array(length, A);
	delete_gate_bootstrapping_ciphertext_array(length, AA);
	delete_gate_bootstrapping_ciphertext_array(length, B);
	delete_gate_bootstrapping_ciphertext_array(length, C);
	delete_gate_bootstrapping_ciphertext_array(length, D);
	delete_gate_bootstrapping_ciphertext_array(length, E);
	delete_gate_bootstrapping_ciphertext_array(length, F);
}



void newRealDiv(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
	
	LweSample* temp = new_gate_bootstrapping_ciphertext_array(3, bk->params);
	LweSample* A = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* B = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* QR = new_gate_bootstrapping_ciphertext_array(2*length, bk->params);
	LweSample* D = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* C = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* Q = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* DD = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* R = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* E0 = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	
	newABS(A, a, length, bk);
	newABS(B, b, length, bk);

	for(int i = 0; i < length; i++){
		bootsCOPY(&QR[i], &A[i], bk);
		bootsCOPY(&D[i], &B[i], bk);
		bootsCONSTANT(&QR[length + i], 0, bk);
        bootsCONSTANT(&C[i], 0, bk);
    }
	

	HomRShift(QR, QR, 2*length, 1, bk);

	for(int s = 1; s < length; s++){
		HomRShift(QR, QR, 2*length, 1, bk);

		for(int i = 0; i < length; i++){
			bootsCOPY(&R[i], &QR[length+i], bk);}

		newCompS(&temp[0], R, D, length, bk);
		bootsNOT(&temp[1], &temp[0], bk);
		bootsCOPY(&QR[0], &temp[1], bk);

		for(int i = 0; i < length; i++){
			bootsAND(&DD[i], &D[i], &temp[1], bk);}

		newSUB(R, R, DD, length, bk);

		for(int i = 0; i < length; i++){
			bootsCOPY(&QR[length+i], &R[i], bk);}
		

	}
		


	for(int s = length; s < length*3/2; s++){
		HomRShift(QR, QR, 2*length, 1, bk);
		HomRShift(R, R, length, 1, bk);

		newCompS(&temp[0], R, D, length, bk);
		bootsNOT(&temp[1], &temp[0], bk);
		bootsCOPY(&QR[0], &temp[1], bk);

		for(int i = 0; i < length; i++){
			bootsAND(&DD[i], &D[i], &temp[1], bk);}

		newSUB(R, R, DD, length, bk);
	}


	for(int i = 0; i < length; i++){
		bootsCOPY(&Q[i], &QR[i], bk);}
	
	bootsXOR(&C[0], &a[length-1], &b[length-1], bk);

	for(int i = 0; i < length; i++){
		bootsXOR(&Q[i], &Q[i], &C[0], bk);
    }

	newADD(res, Q, C, length, bk);
	
	
	
	
	delete_gate_bootstrapping_ciphertext_array(3, temp);
	delete_gate_bootstrapping_ciphertext_array(length, A);
	delete_gate_bootstrapping_ciphertext_array(length, B);
	delete_gate_bootstrapping_ciphertext_array(2*length, QR);
	delete_gate_bootstrapping_ciphertext_array(length, D);
	delete_gate_bootstrapping_ciphertext_array(length, DD);
	delete_gate_bootstrapping_ciphertext_array(length, Q);
	delete_gate_bootstrapping_ciphertext_array(length, R);
	delete_gate_bootstrapping_ciphertext_array(length, C);
	delete_gate_bootstrapping_ciphertext_array(length, E0);
}

void newMult(LweSample* res, const LweSample* a, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    LweSample* tmp = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* A = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* B = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* sign = new_gate_bootstrapping_ciphertext(bk->params);

	newABS(A, a, length, bk);
	newABS(B, b, length, bk);
    for(int i = 0; i < length; i++){
        bootsAND(&res[i], &A[i], &B[0], bk);
    }
    

    for(int i = 1; i < length; i++){
        for (int j = 0; j < length-i; j++){
            bootsAND(&tmp[j], &A[j], &B[i], bk);
        }
        newADD(res+i, res+i, tmp, length-i, bk);
    }

    // bootsCONSTANT(&res[length-1], 0, bk);

	bootsXOR(sign, &a[length-1], &b[length-1], bk);


    for(int i = 0; i < length; i++){
        bootsXOR(&res[i], &res[i], sign, bk);
        bootsCONSTANT(&tmp[i], 0, bk);
    }
    bootsCOPY(&tmp[0], sign, bk);
    newADD(res, res, tmp, length, bk);

}



float cov(float *plain_a, float *plain_b, int32_t N) {

    float A[N], B[N];

    for (int i = 0; i < N; i++) { 
        A[i] = plain_a[i];
        B[i] = plain_b[i];
    }
  
    float ab_sum = 0;
	float a_sum = 0;
	float b_sum = 0;
	
    float ab[N];
	for (int i = 0; i < N; i++) {
		ab[i] = (A[i]*B[i]);
		ab_sum += ab[i];
	}

	for (int i = 0; i < N; i++) {
		a_sum += A[i];
		b_sum += B[i];
	}

	float mult_asum_bsum = a_sum*b_sum;
    float temp1 = mult_asum_bsum / N;
    float temp2 = ab_sum - temp1;	
    float res_ans = temp2 / N;

    return res_ans;
    
}


// void BootsSort(LweSample* res1, LweSample* res2, LweSample* a, LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
// 	LweSample* temp = new_gate_bootstrapping_ciphertext_array(2, bk->params);

// 	HomCompLE(&temp[0], a, b, length, bk);
// 	bootsNOT(&temp[1], &temp[0], bk);

// 	for (int i =0; i < length; i++) {
// 		bootsMUX(&res1[i], &temp[0], &a[i], &b[i], bk);
// 		bootsMUX(&res2[i], &temp[1], &a[i], &b[i], bk);
// 	}

// 	delete_gate_bootstrapping_ciphertext_array(2, temp);	
// }


// biggest one is the first one
// matrix is nx1 matrix
void newMaxValue(LweSample* res, vector<vector<LweSample*>> matrix, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    LweSample* temp = new_gate_bootstrapping_ciphertext_array(2, bk->params);
	LweSample* result1 = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	LweSample* result2 = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    for (int h = 0; h < matrix.size() - 1; h++) {
        for (int i = 0; i < matrix.size() - 1 - h; i++) {
            newCompLE(&temp[0], matrix[i][0], matrix[i+1][0], length, bk); // matrix[i] >= matrix[i+1] -> temp[0] = 1
            bootsNOT(&temp[1], &temp[0], bk);
            for(int k = 0; k < length; k++){
                bootsMUX(&result1[k], &temp[0], &matrix[i][0][k], &matrix[i+1][0][k], bk); // max
                bootsMUX(&result2[k], &temp[1], &matrix[i][0][k], &matrix[i+1][0][k], bk); // min
                bootsCOPY(&matrix[i][0][k], &result1[k], bk);
                bootsCOPY(&matrix[i+1][0][k], &result2[k], bk);
            }
        }
    }

    // W[][] = [v₁, v₂, ..., vr]   // reduced basis
    for(int k = 0; k < length; k++){
        bootsCOPY(&res[k], &matrix[0][0][k], bk);
    }

    // cout << "maxvalue is first number" << endl;
    // decrypting_code(matrix, length);

	delete_gate_bootstrapping_ciphertext_array(2, temp);
	delete_gate_bootstrapping_ciphertext_array(length, result1);
    delete_gate_bootstrapping_ciphertext_array(length, result2);
}


// // biggest one is the latest one
// void newMaxValue(LweSample* res, vector<vector<LweSample*>> matrix, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
// 	LweSample* result1 = new_gate_bootstrapping_ciphertext_array(length, bk->params);
// 	LweSample* result2 = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	
//     for (int i = 0; i < matrix.size(); i++) {
// 	    for (int j = 0; j < matrix[0].size()-1; j++) {
//             // HomMax(result, ciphertext1[j], ciphertext1[j+1], length, bk);
//             // HomMin(result2, ciphertext1[j], ciphertext1[j+1], length, bk);
// 			BootsSort(result1, result2, matrix[j][0], matrix[j+1][0], length, bk);

//             for(int k = 0; k < length; k++){
//                 bootsCOPY(&matrix[j][0][k], &result2[k], bk);
//                 bootsCOPY(&matrix[j+1][0][k], &result1[k], bk);
//             }
//         }
//     }

//     cout << "maxvalue is last number" << endl;
//     decrypting_code(matrix, length);

//     for (int i = 0; i < length; i++)
//         bootsCOPY(&res[i], &matrix[matrix.size() -1][0][i], bk);

// 	delete_gate_bootstrapping_ciphertext_array(length, result1);
//     delete_gate_bootstrapping_ciphertext_array(length, result2);
// }	
