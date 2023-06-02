#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <iostream>
#include <math.h>
#include "add.h"
#include <tfhe/tfhe_gate_bootstrapping_functions.h>
//#include "newmatrix.cpp"
//#include "newmatrix.h"

////////////////////////////////////////
/*vector<vector<float>> res is matrix*/

using namespace std;

int main() {
    printf("============ Step3 start ============\n");
    //reads the secret key from file
    FILE* secret_key = fopen("secret.key","rb");
    TFheGateBootstrappingSecretKeySet* key = new_tfheGateBootstrappingSecretKeySet_fromFile(secret_key);
    fclose(secret_key);

    //if necessary, the params are inside the key
    const TFheGateBootstrappingParameterSet* params = key->params;

    //parameter
    const int length = 32;
    int res_length = 32;
    // const int d = 2;
    // const int r = 1;
    const int dimension = 2;                // at least 2

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


    int row1 = aMatrix.size();
    int col1 = aMatrix[0].size();
    int row2 = bMatrix.size();
    int col2 = bMatrix[0].size();
    int row3 = cMatrix.size();
    int col3 = cMatrix[0].size();
    int row4 = dMatrix.size();
    int col4 = dMatrix[0].size();
    int row5 = 1;
    int col5 = centering.size();


    ////////////////////////////////////////
    /*vector<vector<float>> res is matrix*/

    // read result ciphertext
    ////////////////res matrix
    vector<vector<LweSample*>> answer(row1, vector<LweSample*>(dimension));            ////////////////////        answer      //////////////////////
    for(int i = 0; i < answer.size(); i++) {
        for(int j = 0; j < answer[0].size(); j++) {
            answer[i][j] = new_gate_bootstrapping_ciphertext_array(length, params);
        }
    }

    //import ciphertexts from the answer file
    FILE* answer_data = fopen("answer.data","rb");

    for (int i = 0; i < answer.size(); i++){
        for (int j = 0; j < answer[0].size(); j++){
            for (int k = 0; k < length; k++){
                import_gate_bootstrapping_ciphertext_fromFile(answer_data, &answer[i][j][k], params);
            }
        }
    }

    fclose(answer_data);

    //decrypt and rebuild plaintext answer
    vector<vector<int32_t>> ai(answer.size(),vector<int32_t>(answer[0].size()));
    vector<vector<int32_t>> int_answer(answer.size(),vector<int32_t>(answer[0].size()));
    vector<vector<float>> float_answer(answer.size(),vector<float>(answer[0].size()));

    for (int i = 0; i < answer.size(); i++){
        for (int j = 0; j < answer[0].size(); j++){
            for(int k = 0; k < length; k++){
                ai[i][j] = bootsSymDecrypt(&answer[i][j][k], key);
                int_answer[i][j] |= (ai[i][j]<<k);
            }
        }
    }
    for (int i = 0; i < answer.size(); i++){
        for (int j = 0; j < answer[0].size(); j++){
            float_answer[i][j] = ((float)int_answer[i][j])/((float)pow(2,length/2));
        }
    }

    cout << "### Output Data ###" << endl;

    cout << "=== (tfhe) Matrix _ PCA  === \n"; 

    for (int row = 0; row < answer.size(); row++) {
        for (int col = 0; col < answer[0].size(); col++) {
            cout << float_answer[row][col] << "  ";
        }
        cout << "\n";
    }

    //clean up all pointers
    for (int i = 0; i < answer.size(); i++){
        for (int j = 0; j < answer[0].size(); j++){
            delete_gate_bootstrapping_ciphertext_array(length, answer[i][j]);
        }
    }
    delete_gate_bootstrapping_secret_keyset(key);
    return 0;
}








