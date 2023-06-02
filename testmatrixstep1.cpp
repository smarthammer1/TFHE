#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <iostream>
#include <math.h>
#include "add.h"
#include <tfhe/tfhe_gate_bootstrapping_functions.h>
#include "newfile.cpp"
#include "newfilepca.cpp"
// #include "newmatrix.cpp"
// #include "newmatrix.h"


using namespace std;


int main() {
    printf("============ Step1 start ============\n");
    //generate a keyset
    const int minimum_lambda = 110;
    TFheGateBootstrappingParameterSet* params = new_default_gate_bootstrapping_parameters(minimum_lambda);

    //generate a random key
    uint32_t seed[] = { 314, 1592, 657 };
    tfhe_random_generator_setSeed(seed,3);
    TFheGateBootstrappingSecretKeySet* key = new_random_gate_bootstrapping_secret_keyset(params);

    //export the secret key to file for later use
    FILE* secret_key = fopen("secret.key","wb");
    export_tfheGateBootstrappingSecretKeySet_toFile(secret_key, key);
    fclose(secret_key);

    //export the cloud key to a file (for the cloud)
    FILE* cloud_key = fopen("cloud.key","wb");
    export_tfheGateBootstrappingCloudKeySet_toFile(cloud_key, &key->cloud);
    fclose(cloud_key);

    //you can put additional instructions here!!
    //////////////////////////////////////////////////////////////
    const int length = 32;
    // const int d = 2;
    // const int r = 1;
    const int number_of_iteration = 3;      // at least 3
    float num = 3.6;
    const int dimension = 2;                // at least 2

    ////generate ciphertext from plaintext
    vector<vector<float>> aMatrix = {{1, 13, 4},{-3, -39, -12}, {1.6, 20.8, 6.4}};
    // vector<vector<float>> aMatrix = {{1.3, 0, 7.5}, {1, -3.6, -4}};
    // vector<vector<float>> aMatrix = {{2.5, 0}, {4.1, -7}};
    // vector<vector<float>> aMatrix = {{1.5, 1.3}};
    // vector<vector<float>> bMatrix = {{1.3},{3}};

// /*(2x3)*/    vector<vector<float>> bMatrix = {{2.3, -4.5, 7.3}, {0, 1, -2}};
// /*(4x3)*/    vector<vector<float>> bMatrix = {{1.5, 2.7, 0}, {3, -4, 5.4}, {2.3, 1, -4.5}, {6.9, -5.5, 2}};
// /*(4x4)*/    vector<vector<float>> bMatrix = {{1.5, 2.7, 0, 5}, {3, -4, 5.4, -8}, {2.3, 1, -4.5, 9}, {6.9, -5.5, 2, 7}};
// /*(6x5)*/    vector<vector<float>> bMatrix = {{1.5, 2.7, 0, 5, -7.6}, {3, -4, 0, 3.6, 5.4}, {2.3, 1, 6, 8.2, -4.5}, {6.9, -5.5, 0, -6.2, 2}, {1.5, 0.3, 2.1, -0.9, 4.3}, {1.2, 3.4, 2.5, 0.7, 4.1}};
// /*(6x6)*/    vector<vector<float>> bMatrix = {{1.5, 2.7, 0, 5, -7.6, 8}, {3, -4, 0, 3.6, 5.4, -8}, {2.3, 1, 6, 8.2, -4.5, 9}, {6.9, -5.5, 0, -6.2, 2, 7}, {1.5, 0.3, 2.1, -0.9, 4.3, -1.7}, {1.2, 3.4, 2.5, 0.7, 4.1, 1.9}};
// /*(8x8)*/    vector<vector<float>> bMatrix = {{1.5, 2.7, 0, 5, -7.6, 8, 2.2, -0.5}, {3, -4, 0, 3.6, 5.4, -8, -1.2, 0.6}, {2.3, 1, 6, -0.9, 4.3, 8.2, -4.5, 9}, {6.9, -5.5, 0, -6.2, 4.6, -2.8, 2, 7}, {1.5, 0.3, 2.1, -0.9, 4.1, 1.9, 4.3, -1.7}, {1.2, 3.4, 2.5, -0.4, 2.7, 0.7, 4.1, 1.9}, {7.2, -1.9, 0.6, 2.4, 1.5, 2.7, 3.2, 4.1}, {-1.6, -3.9, 2.7, 0, 2.9, 5.6, 6.8, -1.1}};       
// /*(2x10)*/    vector<vector<float>> bMatrix = {{1.5, 2.7, 0, 5, -7.6, 8, 3, -4, 0, -8}, {2.3, 1, 6, 8.2, -4.5, 9, 6.9, 0, -6.2, 2}};
/*(10x10)*/    vector<vector<float>> bMatrix = {{2.5, 7.6, 2.3, -1.3, 0.7, -4.8, 3.1, 0, -6, -0.5}, {3.8, 0.2, -7.4, 5.6, -2.9, 6, 4, 8, -1.3, 8.5}, {4.5, 7.8, -2.3, -6.7, 1.2, 5.6, 3.4, -9, 0.1, 8}, {4.3, -7.5, 9.8, 2.5, -5.7, 8.9, 6.3, -3.1, 1.5, 0.6}, {-1, 3.2, 6.6, -8.8, 10, 2.2, 5.5, -7.7, 0, 4.4}, {8.7, -6.3, 4.9, 1.3, -5, 2.7, 0.9, 9.2, -7, 3}, {3.8, 0, 7.4, -5.6, 2.9, 9, 6.8, -4.7, 8.2, 1.3}, {2.4, 1.2, -3.6, 5.8, 4.1, 9.7, 6.5, -8.9, 7.3, -0.9}, {1.5, 2.7, 0, 5, -7.6, 8, 3, -4, 0, -8}, {2.3, 1, 6, 8.2, -4.5, 9, 6.9, 0, -6.2, 2}};

    vector<vector<float>> cMatrix = {{3.1, 1}, {-5, 3.8}};
    // vector<vector<float>> cMatrix = {{3.2}, {-1}, {0.4}};
    vector<vector<float>> dMatrix = {{2.3, -4.5}, {7.3, 0}, {1, -2}};
    vector<float> centering = {{3.3}, {-3.5}};
    float number = 3.9;


    // encode float to plaintext
    vector<vector<int32_t>> encode_plain1(aMatrix.size(), vector<int32_t>(aMatrix[0].size()));
    for (int i = 0; i < aMatrix.size(); i++){
        for (int j = 0; j < aMatrix[0].size(); j++){
            encode_plain1[i][j] = (int32_t)(aMatrix[i][j]*pow(2,length/2));
        }
    }
    vector<vector<int32_t>> encode_plain2(bMatrix.size(), vector<int32_t>(bMatrix[0].size()));
    for (int i = 0; i < bMatrix.size(); i++){
        for (int j = 0; j < bMatrix[0].size(); j++){
            encode_plain2[i][j] = (int32_t)(bMatrix[i][j]*pow(2,length/2));
        }
    }
    vector<vector<int32_t>> encode_plain3(cMatrix.size(), vector<int32_t>(cMatrix[0].size()));
    for (int i = 0; i < cMatrix.size(); i++){
        for (int j = 0; j < cMatrix[0].size(); j++){
            encode_plain3[i][j] = (int32_t)(cMatrix[i][j]*pow(2, length/2));
        }
    }
    vector<vector<int32_t>> encode_plain4(dMatrix.size(), vector<int32_t>(dMatrix[0].size()));
    for (int i = 0; i < dMatrix.size(); i++){
        for (int j = 0; j < dMatrix[0].size(); j++){
            encode_plain4[i][j] = (int32_t)(dMatrix[i][j]*pow(2, length/2));
        }
    }
    vector<int32_t> encode_plain5(centering.size());
    for (int i = 0; i < centering.size(); i++){
        encode_plain5[i] = (int32_t)(centering[i]*pow(2, length/2));
    }
    int32_t encode_plain6 = (int32_t)(number*pow(2, length/2));


    // plaintext1 to ciphertext1: plaintext bit-by-bit encryption to ciphertext
    int row1 = aMatrix.size();
    int col1 = aMatrix[0].size();
    vector<vector<LweSample*>> cipher1(row1,vector<LweSample*>(col1));
    for (int i = 0; i < aMatrix.size(); i++){
        for (int j = 0; j < aMatrix[0].size(); j++){
            cipher1[i][j] = new_gate_bootstrapping_ciphertext_array(length, params);
        }
    }

    int row2 = bMatrix.size();
    int col2 = bMatrix[0].size();
    vector<vector<LweSample*>> cipher2(row2,vector<LweSample*>(col2));
    for (int i = 0; i < bMatrix.size(); i++){
        for (int j = 0; j < bMatrix[0].size(); j++){
            cipher2[i][j] = new_gate_bootstrapping_ciphertext_array(length, params);
        }
    }
    
    int row3 = cMatrix.size();
    int col3 = cMatrix[0].size();
    vector<vector<LweSample*>> cipher3(row3, vector<LweSample*>(col3));
    for (int i = 0; i < cMatrix.size(); i++){
        for (int j = 0; j < cMatrix[0].size(); j++){
            cipher3[i][j] = new_gate_bootstrapping_ciphertext_array(length, params);
        }
    }

    int row4 = dMatrix.size();
    int col4 = dMatrix[0].size();
    vector<vector<LweSample*>> cipher4(row4, vector<LweSample*>(col4));
    for (int i = 0; i < dMatrix.size(); i++){
        for (int j = 0; j < dMatrix[0].size(); j++){
            cipher4[i][j] = new_gate_bootstrapping_ciphertext_array(length, params);
        }
    }

    int row5 = 1;
    int col5 = centering.size();
    vector<LweSample*> cipher5(col5);
    for (int i = 0; i < centering.size(); i++){
        cipher5[i] = new_gate_bootstrapping_ciphertext_array(length, params);
    }

    LweSample* cipher6 = new_gate_bootstrapping_ciphertext_array(length, params);


    for (int i = 0; i < aMatrix.size(); i++){
        for (int j = 0; j < aMatrix[0].size(); j++){
            for (int k = 0; k < length; k++){
                bootsSymEncrypt(&cipher1[i][j][k], (encode_plain1[i][j]>>k)&1, key);
            }
        }
    }
    for (int i = 0; i < bMatrix.size(); i++){
        for (int j = 0; j < bMatrix[0].size(); j++){
            for (int k = 0; k < length; k++){
                bootsSymEncrypt(&cipher2[i][j][k], (encode_plain2[i][j]>>k)&1, key);
            }
        }
    }
    for (int i = 0; i < cMatrix.size(); i++){
        for (int j = 0; j < cMatrix[0].size(); j++){
            for (int k = 0; k < length; k++){
                bootsSymEncrypt(&cipher3[i][j][k], (encode_plain3[i][j]>>k)&1, key);
            }
        }
    }
    for (int i = 0; i < dMatrix.size(); i++){
        for (int j = 0; j < dMatrix[0].size(); j++){
            for (int k = 0; k < length; k++){
                bootsSymEncrypt(&cipher4[i][j][k], (encode_plain4[i][j]>>k)&1, key);
            }
        }
    }
    for (int i = 0; i < centering.size(); i++){
        for (int k = 0; k < length; k++){
            bootsSymEncrypt(&cipher5[i][k], (encode_plain5[i]>>k)&1, key);
        }
    }
    for (int k = 0; k < length; k++){
        bootsSymEncrypt(&cipher6[k], (encode_plain6>>k)&1, key);
    }


    // /*
    // print before result_decrypting encoded plaintext
    // int32_t tmp;
    // for (int k = 0; k < length; k++){
    //     // tmp = bootsSymDecrypt(&cipher1[1][0][k], key);
    //     cout <<  ((encode_plain1[0][0]>>k)&1) << "\n";
    // }
    // */


    //show algorithm type
    cout << "########### Matrix Algorithm ###########" << endl;

    //show your data
    cout << "### Input Data ###" << endl;

    printf("%d-by-%d matrix\n", row2, col2);

    printf("=== Matrix a ===\n");
    printMatrix(aMatrix);
    printf("\n");

    // printf("=== Matrix b ===\n");
    // printMatrix(bMatrix);
    // printf("\n");    
    
    // printf("=== Matrix c ===\n");
    // printMatrix(cMatrix);
    // printf("\n");

    // printf("=== Matrix d ===\n");
    // printMatrix(dMatrix);
    // printf("\n");    

    // printf("=== Matrix_centering ===\n");
    // for(int i = 0; i < centering.size(); i++){
    //     cout << centering[i] << "  \n";
    // }

    // printf("=== number e ===\n");
    // printMatrix(number);
    // printf("\n";)

    // printf("=== number n ===\n");
    // printf("%f\n",num);
 
    printf("===iteration number===\n");
    printf("%d\n",number_of_iteration);

    printf("=== dimension ===\n");
    printf("%d\n",dimension);


    //export the ciphertexts to a file (for the cloud)
    FILE* cloud_data = fopen("cloud.data","wb");
    for (int i = 0; i < aMatrix.size(); i++){
        for (int j = 0; j < aMatrix[0].size(); j++){
            for (int k = 0; k < length; k++){
                 export_gate_bootstrapping_ciphertext_toFile(cloud_data, &cipher1[i][j][k], params);
            }
        }
    }
    for (int i = 0; i < bMatrix.size(); i++){
        for (int j = 0; j < bMatrix[0].size(); j++){
            for (int k = 0; k < length; k++){
                 export_gate_bootstrapping_ciphertext_toFile(cloud_data, &cipher2[i][j][k], params);
            }
        }
    }
    for (int i = 0; i < cMatrix.size(); i++){
        for (int j = 0; j < cMatrix[0].size(); j++){
            for (int k = 0; k < length; k++){
                export_gate_bootstrapping_ciphertext_toFile(cloud_data, &cipher3[i][j][k], params);
            }
        }
    }
    for (int i = 0; i < dMatrix.size(); i++){
        for (int j = 0; j < dMatrix[0].size(); j++){
            for (int k = 0; k < length; k++){
                export_gate_bootstrapping_ciphertext_toFile(cloud_data, &cipher4[i][j][k], params);
            }
        }
    }
    for (int i = 0; i < centering.size(); i++){
        for (int k = 0; k < length; k++){
            export_gate_bootstrapping_ciphertext_toFile(cloud_data, &cipher5[i][k], params);
        }
    }
    for (int k = 0; k < length; k++){
        export_gate_bootstrapping_ciphertext_toFile(cloud_data, &cipher6[k], params);
    }


    fclose(cloud_data);


    //clean up all pointers
    for (int i = 0; i < aMatrix.size(); i++){
        for (int j = 0; j < aMatrix[0].size(); j++){
            delete_gate_bootstrapping_ciphertext_array(length, cipher1[i][j]);
        }
    }
    for (int i = 0; i < bMatrix.size(); i++){
        for (int j = 0; j < bMatrix[0].size(); j++){
            delete_gate_bootstrapping_ciphertext_array(length, cipher2[i][j]);
        }
    }
    for (int i = 0; i < cMatrix.size(); i++){
        for (int j = 0; j < cMatrix[0].size(); j++){
            delete_gate_bootstrapping_ciphertext_array(length, cipher3[i][j]);
        }
    }
    for (int i = 0; i < dMatrix.size(); i++){
        for (int j = 0; j < dMatrix[0].size(); j++){
            delete_gate_bootstrapping_ciphertext_array(length, cipher4[i][j]);
        }
    }
    for (int i = 0; i < centering.size(); i++){
        delete_gate_bootstrapping_ciphertext_array(length, cipher5[i]);
    }
    delete_gate_bootstrapping_ciphertext_array(length, cipher6);


    //////////////////////////////////////////////////////////////
    //clean up all pointers
    delete_gate_bootstrapping_secret_keyset(key);
    delete_gate_bootstrapping_parameters(params);
}

