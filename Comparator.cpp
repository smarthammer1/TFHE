#include "add.h"
#include "add.cpp"
#include "Comparator.h"
#include <tfhe/tfhe_gate_bootstrapping_functions.h>
#include <cstdint>
#include <math.h>

#define CLOCKS_PER_MS 1000
extern float g_time;

int64_t encode_to (int64_t res, const int length) {
    switch(length) {
        case 8:
            res = (int8_t) res;
            break;

        case 16:
            res = (int16_t) res;
            break;

        case 32:
            res = (int32_t) res;
            break;

        case 64:
            res = (int64_t) res;
            break;
        
        default:
            cout << "Nothing happend! \n";
            break;
            
    }
return res;
}

int32_t getlength() {

    int32_t len;

    cout << "Enter bit representation length n (8/16/32/64): ";
    cin >> len;

    if (cin.fail()){
        cout << "Wrong Number. Retry! \n"; 
        cin.clear(); 
        cin.ignore(32, '\n'); 
        return getlength(); 
    }
    if (!((len == 8) || (len == 16) || (len == 32) || (len == 64))){
        cout << "Sorry! We provide only n = 8, 16, 32, 64! \n";
        return getlength();
    }
    return len;
}

float getNumber() {

    float retNum;

    cout << "Enter Number: ";
    cin >> retNum;

    if (cin.fail()){
        cout << "Wrong Number. Retry! \n"; 
        cin.clear(); 
        cin.ignore(32, '\n'); 
        return getNumber(); 
    }
    if (retNum < INT32_MIN && retNum > INT32_MAX ){
        cout << "Wrong Number. We recommend your number should be inside " << INT32_MIN << "<= x <= " << INT32_MAX << "\n";
        cin.clear();
        cin.ignore(32, '\n');
        return getNumber();
    }
    return retNum;
}

int32_t getsecparam() {

    int32_t secparam;

    cout << "Enter Security Parameter! (80, 127) : ";
    cin >> secparam;

    if (cin.fail()){
        cout << "Wrong Number. Retry! \n"; 
        cin.clear(); 
        cin.ignore(32, '\n'); 
        return getsecparam(); 
    }
    if (secparam > 127) {
        cout << "TFHE only supports parameters under 127 bit of security! \n";
        return getsecparam();
    }
    else if (secparam > 80 and secparam <= 127) {
        cout << "Security level : 127 bit! \n";
        secparam = 127;
    }
    else if (secparam > 0 and secparam <= 80) {
        cout << "Security level : 80 bit! \n";
        secparam = 80;
    }
    else  {
        cout << "Wrong security parameter! Try again :) \n";
        return getsecparam();
    }

    return secparam;
}

void get_time (float time) {
    g_time = time;
} 

int32_t get_N() {

    float N;

    cout << "Enter Number of Values! : ";
    cin >> N;

    if (cin.fail()){
        cout << "Wrong Number. Retry! \n"; 
        cin.clear(); 
        return get_N(); 
    }
    if (N < 1 || N > INT8_MAX ){
        cout << "Wrong Number. We recommend your number should be inside " << 1 << "<= x <= " << INT8_MAX << "\n";
        cin.clear();
        return get_N();
    }
    if (N != (int)N){
        cout << "Wrong Number. We recommend your number should be an integer inside " << 1 << "<= x <= " << INT8_MAX << "\n";
        cin.clear();
        return get_N();
    }


    return N;
}



void comp_ft(float *plain_a, float *plain_b, LweSample** a, LweSample** b, const int length, int32_t N, int32_t sec_param, const TFheGateBootstrappingCloudKeySet* bk, TFheGateBootstrappingSecretKeySet* key) {
    int sel;
    cout.precision(6);
    
    float plain_a_temp[N], plain_b_temp[N];
    for(int i = 0; i < N; i++){
        plain_a_temp[i] = plain_a[i];
        plain_b_temp[i] = plain_b[i];
    }
    
    LweSample* a_temp[N];
    LweSample* b_temp[N];
    for(int i = 0; i < N; i++){
        a_temp[i]= new_gate_bootstrapping_ciphertext_array(length, bk->params);
        b_temp[i]= new_gate_bootstrapping_ciphertext_array(length, bk->params);
    }

    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < length; j++){
            bootsCOPY(&a_temp[i][j], &a[i][j], bk);
            bootsCOPY(&b_temp[i][j], &b[i][j], bk);
        }
    }


    while(1) {
        cout << "-------------------------------------------------------------------" << "\n";
        cout << "|           0. Gate Bootstrapping Parameters & Runtime            |" << "\n";
        cout << "|                           1. Covariance                         |" << "\n";
        cout << "|                           2. Multiplication                     |" << "\n";
        cout << "|                            -1. Exit                             |" << "\n";
        cout << "-------------------------------------------------------------------" << "\n\n";

        cout << " Please Select the operation number you want to proceed (0 ~ 2, -1) : ";
        cin >> sel;
        if (cin.fail()){
            cout << "Wrong Number. Retry! \n"; 
            cin.clear(); 
            cin.ignore(32, '\n'); 
        }
        
        switch(sel) {
            case 0:
                boot_eval(sec_param, length, bk, key);
                break;
            
            case 1:
                comp_covar(plain_a_temp, plain_b_temp, a_temp, b_temp, length, N, g_time, bk, key);
                break;

            case 2:
                comp_RealMult(plain_a_temp, plain_b_temp, a_temp, b_temp, length, N, g_time, bk, key);
                break;

            case -1:
                exit(1);
                
            default:
                cout << "\n Please select a number inside 0 ~ 2, -1 ! \n";
                break;
        }
    }
    
}

void boot_eval(int sec_param, const int length, const TFheGateBootstrappingCloudKeySet* bk, TFheGateBootstrappingSecretKeySet* key){
    LweSample* tmp1 = new_gate_bootstrapping_ciphertext(bk->params);
    LweSample* tmp2 = new_gate_bootstrapping_ciphertext(bk->params);
    bootsSymEncrypt(tmp1, 0, key);
    bootsSymEncrypt(tmp2, 1, key);
    int mes = 30;
    // LweSample* tmp2 = new_gate_bootstrapping_ciphertext(bk->params);
    LweSample* res = new_gate_bootstrapping_ciphertext(bk->params);
    Torus32 mu = modSwitchToTorus32(1, 8);

    cout << "-------------------------------------------------------------------" << "\n";
    cout << "|           0. Gate Bootstrapping Parameters & Runtime            |" << "\n";
    cout << "-------------------------------------------------------------------" << "\n\n";

    cout << "            Security Parameter :  " << sec_param << "-bits of Security                " << "\n";
    cout << "                   Ciphertext Length :  " << length << "             " << "\n";
    cout << "\n";
    cout << "                         --LWE PARAMS--                                              " << "\n";
    cout << "                             n :  " << bk->params->in_out_params->n << "              " << "\n";
    cout << "                       LWE std :  " << bk->params->in_out_params->alpha_min << "      " << "\n";
    cout << "\n";
    cout << "                        --TLWE PARAMS--                                              " << "\n";
    cout << "                             N :  " << bk->params->tgsw_params->tlwe_params->N << "   " << "\n";
    cout << "                             k :  " << bk->params->tgsw_params->tlwe_params->k << "   " << "\n";
    cout << "                      TLWE std :  " << bk->params->tgsw_params->tlwe_params->alpha_min << "   " << "\n";
    cout << "\n";
    cout << "                        --TGSW PARAMS--                                              " << "\n";
    cout << "                            Bg :  " << bk->params->tgsw_params->Bg << "               " << "\n";
    cout << "                             l :  " << bk->params->tgsw_params->l << "                " << "\n";
    cout << "\n";
    float time = -clock();
    for(int i = 0; i < mes; i++){
        bootsAND(res, tmp1, tmp2, bk);
    }
    time += clock();
    time = time/(mes * CLOCKS_PER_MS);
    cout << "        Time per Bootstrapping :  " << time << " ms \n";
    get_time(time);


    delete_gate_bootstrapping_ciphertext(tmp1);
    delete_gate_bootstrapping_ciphertext(tmp2);
    delete_gate_bootstrapping_ciphertext(res);
}

void comp_covar(float *plain_a, float *plain_b, LweSample** a, LweSample** b, const int length, int32_t N, float bs_time, const TFheGateBootstrappingCloudKeySet* bk, TFheGateBootstrappingSecretKeySet* key) {
    LweSample* res1 = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* res2 = new_gate_bootstrapping_ciphertext_array(length, bk->params);


    float plain_a_temp[N], plain_b_temp[N];
    for(int i = 0; i < N; i++){
        plain_a_temp[i] = plain_a[i];
        plain_b_temp[i] = plain_b[i];
    }

    LweSample* a_temp[N];
    LweSample* b_temp[N];

    for(int i = 0; i < N; i++){
        a_temp[i]= new_gate_bootstrapping_ciphertext_array(length, bk->params);
        b_temp[i]= new_gate_bootstrapping_ciphertext_array(length, bk->params);
    }

    for(int i = 0; i < N; i++){
        for(int j = 0; j < length; j++){
            bootsCOPY(&a_temp[i][j], &a[i][j], bk);
            bootsCOPY(&b_temp[i][j], &b[i][j], bk);
        }
    }


    cout << "\n";
    cout << "-------------------------------------------------------------------" << "\n";
    cout << "|                           1. Covariance                         |" << "\n";
    cout << "-------------------------------------------------------------------" << "\n";
    cout << "\n                    plaintext 1  :  "; for(const auto &value: plain_a_temp){std::cout << value << ' ';}  
    cout << "\n                    plaintext 2  :  "; for(const auto &value: plain_b_temp){std::cout << value << ' ';} 
    cout << "\n\n       We will evaluate the time taken for " << length << "-bit Covariance! \n\n";
    
    cout << "                newCovar           |              HomCovar                 " << "\n";
    cout << "                                   |\n";
    cout << "  # of op :    "; printf("%.1fn^2+%.1fn-%d", 1.5*N+12, 13.5*N+25.5, 4*N+4 );  
    cout << "   |           ";   printf("%.1fn^2+%.1fn-%d\n", 3*N+22.5, 21*N+21.5, 14*N+23 );
    cout << "     Est  :   " << bs_time*((float) (1.5*N+12)*pow(length, 2) + (float) (13.5*N+25.5)*length -(4*N+4) ) << " ms            |          " 
         << bs_time*((float) (3*N+22.5)*pow(length, 2) + (float) (21*N+21.5)*length -(14*N+23) ) << " ms \n";
    
    

    float time = -clock();
    newCovar(res1, a_temp, b_temp, length, N, bk);
    time += clock();
    time = time/CLOCKS_PER_MS;


    float time2 = -clock();
    HomCovar(res2, a_temp, b_temp, length, N, bk);
    time2 += clock();
    time2 = time2/CLOCKS_PER_MS;

    cout << "    Real  :   " << time << " ms         |            " << time2 << " ms \n";
    cout << "\n";
    cout << "                     " << time2/time << "x times of Speedup! \n";
    cout << "\n";
    int32_t answer = 0;
    for (int i=0; i<length; i++) {
        int32_t ai = bootsSymDecrypt(&res1[i], key);
        answer |= (ai<<i);
    }
    float float_answer = ((float)answer)/((float)pow(2.0,16.0));

    // answer = bootsSymDecrypt(compres1, key);
    
    int32_t answer2 = 0;
    for (int i=0; i<length; i++) {
        int32_t ai2 = bootsSymDecrypt(&res2[i], key);
        answer2 |= (ai2<<i);
    }
    float float_answer2 = ((float)answer2)/((float)pow(2.0,16.0));
    
    
    printf("                 Plaintext Covariance :  cov( (");
    for(const auto &value: plain_a_temp){std::cout << value << ' ';} 
    printf("), ("); 
    for(const auto &value: plain_b_temp){std::cout << value << ' ';} 
    printf(") ) = ");  
    cout <<  cov(plain_a_temp, plain_b_temp, N) << "\n";
    
    cout << "                 newCovar Covariance :  " << float_answer << " \n";
    cout << "                 HomCovar Covariance :  " << float_answer2 << " \n\n";
    
    delete_gate_bootstrapping_ciphertext_array(length, res1);
    delete_gate_bootstrapping_ciphertext_array(length, res2);
}



void comp_RealMult(float *plain_a, float *plain_b, LweSample** a, LweSample** b, const int length, int32_t N, float bs_time, const TFheGateBootstrappingCloudKeySet* bk, TFheGateBootstrappingSecretKeySet* key) {
    LweSample* res1 = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* res2 = new_gate_bootstrapping_ciphertext_array(length, bk->params);

    float plain_a_temp[N], plain_b_temp[N];
    for(int i = 0; i < N; i++){
        plain_a_temp[i] = plain_a[i];
        plain_b_temp[i] = plain_b[i];
    }

    LweSample* a_temp[N];
    LweSample* b_temp[N];

    for(int i = 0; i < N; i++){
        a_temp[i]= new_gate_bootstrapping_ciphertext_array(length, bk->params);
        b_temp[i]= new_gate_bootstrapping_ciphertext_array(length, bk->params);
    }

    for(int i = 0; i < N; i++){
        for(int j = 0; j < length; j++){
            bootsCOPY(&a_temp[i][j], &a[i][j], bk);
            bootsCOPY(&b_temp[i][j], &b[i][j], bk);
        }
    }


    cout << "\n";
    cout << "-------------------------------------------------------------------" << "\n";
    cout << "|                        2. Multiplication                       |" << "\n";
    cout << "-------------------------------------------------------------------" << "\n";
    cout << "\n                    plaintext 1  :  "; for(const auto &value: plain_a_temp){std::cout << value << ' ';}  
    cout << "\n                    plaintext 2  :  "; for(const auto &value: plain_b_temp){std::cout << value << ' ';} 
    cout << "     We will evaluate the time taken for " << length << "-bit Multiplication! \n\n";
    
    cout << "                newMultiReal          |             HomMultiReal                 " << "\n";
    cout << "                                 |\n";
    cout << "  # of op : 3/2n^2 + 15/2n - 1   |          3n^2 + 6n - 5                 \n";  
    cout << "     Est  :   " << bs_time*((float) 1.5*pow(length, 2) + 7.5*length -1) << " ms         |            " << bs_time*((float) 3*pow(length, 2) + 6*length - 5) << " ms \n";

    float time = -clock();
    
    newMultiReal(res1, a_temp, b_temp, length, N, bk);
    time += clock();
    time = time/CLOCKS_PER_MS;

    float time2 = -clock();
    HomMultiReal(res2, a_temp, b_temp, length, N, bk);
    time2 += clock();
    time2 = time2/CLOCKS_PER_MS;
    
    cout << "    Real  :   " << time << " ms         |            " << time2 << " ms \n";
    cout << "\n";
    cout << "                     " << time2/time << "x times of Speedup! \n";
    cout << "\n";

    int64_t answer = 0;
    for (int i=0; i<length; i++) {
        int64_t ai = bootsSymDecrypt(&res1[i], key);
        answer |= (ai<<i);
    }
    float float_answer = ((float)answer)/((float)pow(2.0,16.0));

    // answer = bootsSymDecrypt(compres1, key);

    int64_t answer2 = 0;
    for (int i=0; i<length; i++) {
        int64_t ai2 = bootsSymDecrypt(&res2[i], key);
        answer2 |= (ai2<<i);
    }
    float float_answer2 = ((float)answer2)/((float)pow(2.0,16.0));

    int64_t res_ans = 0;
    res_ans = encode_to(plain_a*plain_b, length);

    cout << "        Plaintext Multiplication :  " << plain_a << " * "<< plain_b << " = " << res_ans << "                " << "\n";
    cout << "          newMultiReal Multiplication :  " << float_answer << " \n";
    cout << "          HomMultiReal Multiplication :  " << float_answer2 << " \n\n";

    delete_gate_bootstrapping_ciphertext_array(length, res1);
    delete_gate_bootstrapping_ciphertext_array(length, res2);
}