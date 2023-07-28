#include "add.h"
#include "add.cpp"
#include "Comparator.h"
#include <tfhe/tfhe_gate_bootstrapping_functions.h>
#include <cstdint>
#include <math.h>
#include <iostream>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

using namespace std;

void printMatrix(vector<vector<float>> matrix) {
    for (int row = 0; row < matrix.size(); row++) {
        for (int col = 0; col < matrix[0].size(); col++) {
            cout << matrix[row][col] << "  ";            
        }
        cout << "\n";
    }
}

vector<vector<float>> print_matrix(vector<vector<float>> matrix) {
    vector<vector<float>> plainres(matrix.size(), vector<float>(matrix.size()));
    for (int row = 0; row < matrix.size(); row++) {
        for (int col = 0; col < matrix[0].size(); col++) {
            plainres[row][col] = matrix[row][col];            
        }
    }
    return plainres;
}

vector<vector<float>> AddMatrix(vector<vector<float>> aMatrix, vector<vector<float>> bMatrix) {
    vector<vector<float>> plainres(aMatrix.size(), vector<float>(aMatrix[0].size()));
    for (int row = 0; row < aMatrix.size(); row++) {
        for (int col = 0; col < aMatrix[0].size(); col++) {
            plainres[row][col] = aMatrix[row][col] + bMatrix[row][col];
        }
    }
    return plainres;
}

vector<vector<float>> SubtractMatrix(vector<vector<float>> aMatrix, vector<vector<float>> bMatrix) {
    vector<vector<float>> plainres(aMatrix.size(), vector<float>(aMatrix[0].size()));
    for (int row = 0; row < aMatrix.size(); row++) {
        for (int col = 0; col < aMatrix[0].size(); col++) {
            plainres[row][col] = aMatrix[row][col] - bMatrix[row][col];
        }
    }
    return plainres;
}

vector<vector<float>> MultiplyMatrix(vector<vector<float>> aMatrix, vector<vector<float>> bMatrix) {
    vector<vector<float>> plainres(aMatrix.size(), vector<float>(bMatrix[0].size()));
    for (int row = 0; row < aMatrix.size(); row++) {
        for (int col = 0; col < bMatrix[0].size(); col++) {
            float tmp = 0;
            // Multiply the row of A by the column of B to get the row, column of product.
            for (int inner = 0; inner < aMatrix[0].size(); inner++) {
                tmp += aMatrix[row][inner] * bMatrix[inner][col];
            }
            plainres[row][col] = tmp;
        }
    }
    return plainres;
}

vector<vector<float>> TransposeMatrix(vector<vector<float>> matrix) {
    vector<vector<float>> plainres(matrix[0].size(), vector<float>(matrix.size()));
    for (int row = 0; row < matrix[0].size(); row++) {
        for (int col = 0; col < matrix.size(); col++) {
            plainres[row][col] = matrix[col][row];
        }
    }
    return plainres;
}

// void InverseMatrix(vector<vector<float>> matrix, const int number_of_iteration) {
//     // number_of_iteration = 1
//     vector<vector<float>> A(matrix.size(), vector<float>(matrix[0].size()));
//     vector<vector<float>> transposeA(matrix[0].size(), vector<float>(matrix.size()));
//     vector<vector<float>> multMatrix(matrix.size(), vector<float>(matrix.size()));
//     float traceA = 0;
//     vector<vector<float>> identityMatrix(matrix.size(), vector<float>(matrix.size()));
//     vector<vector<float>> temp(matrix.size(), vector<float>(matrix.size()));
//     vector<vector<float>> ab(matrix.size(), vector<float>(matrix.size()));
//     vector<vector<float>> res_ans(matrix[0].size(), vector<float>(matrix.size()));
    
//     for (int i = 0; i < matrix.size(); i++) {
//         for (int j = 0; j < matrix.size(); j++) {
//             multMatrix[i][j] = 0;
//             ab[i][j] = 0;
//             identityMatrix[i][j] = 0;
//         }
//     }
//     for(int i = 0; i < matrix.size(); i++) {
//         identityMatrix[i][i] = 1;
//     }
//     for (int i = 0; i < matrix[0].size(); i++) {
//         for (int j = 0; j < matrix.size(); j++) {
//             res_ans[i][j] = 0;
//         }
//     }

//     for (int i = 0; i < matrix.size(); i++) {
//         for (int j = 0; j < matrix[0].size(); j++) {
//             A[i][j] = matrix[i][j];
//             transposeA[j][i] = matrix[i][j];
//         }
//     }
//     for (int i = 0; i < matrix.size(); i++) {
//         for (int k = 0; k < matrix.size(); k++) {
//             for (int j = 0; j < matrix[0].size(); j++) {
//                 multMatrix[i][k] += (A[i][j])*(transposeA[j][k]);
//             }
//         }
//     }
//     for (int i = 0; i < matrix.size(); i++) {
//         traceA += multMatrix[i][i];
//     }
//     for (int i = 0; i < matrix[0].size(); i++) {
//         for (int j = 0; j < matrix.size(); j++) {
//             transposeA[i][j] = (transposeA[i][j])*(1/traceA);
//         }
//     }

//     for(int n = 0; n < number_of_iteration; n++){
//         for (int i = 0; i < matrix.size(); i++) {
//             for (int k = 0; k < matrix.size(); k++) {
//                 for (int j = 0; j < matrix[0].size(); j++) {
//                     ab[i][k] += (A[i][j])*(transposeA[j][k]);
//                 }
//             }
//         }
//         for (int i = 0; i < matrix.size(); i++) {
//             for (int j = 0; j < matrix.size(); j++) {
//                 temp[i][j] = (identityMatrix[i][j] - ab[i][j]);
//             }
//         }
//         for (int i = 0; i < matrix[0].size(); i++) {
//             for (int k = 0; k < matrix.size(); k++) {
//                 for (int j = 0; j < matrix.size(); j++) {
//                     res_ans[i][k] += (transposeA[i][j])*(temp[j][k]);
//                 }
//             }
//         }
//         for (int i = 0; i < matrix[0].size(); i++){
//             for (int k = 0; k < matrix.size(); k++) {
//                 transposeA[i][k] = res_ans[i][k];
//             }
//         }
//     }

//     for (int i = 0; i < matrix[0].size(); i++) {
//         for (int j = 0; j < matrix.size(); j++) {
//             cout << transposeA[i][j] << " ";
//         }
//         cout << "\n";
//     }
// }





void decrypting_code(vector<vector<LweSample*>> res, const int length) {
    //reads the secret key from file
    FILE* secret_key = fopen("secret.key","rb");
    TFheGateBootstrappingSecretKeySet* key = new_tfheGateBootstrappingSecretKeySet_fromFile(secret_key);
    fclose(secret_key);

    //if necessary, the params are inside the key
    const TFheGateBootstrappingParameterSet* params = key->params;

    //decrypt and rebuild plaintext answer
    vector<vector<int32_t>> ai(res.size(),vector<int32_t>(res[0].size()));
    vector<vector<int32_t>> int_res(res.size(),vector<int32_t>(res[0].size()));
    vector<vector<float>> float_res(res.size(),vector<float>(res[0].size()));
    for (int i = 0; i < res.size(); i++){
        for (int j = 0; j < res[0].size(); j++){
            for(int k = 0; k < length; k++){
                ai[i][j] = bootsSymDecrypt(&res[i][j][k], key);
                int_res[i][j] |= (ai[i][j]<<k);
            }
        }
    }
    for (int i = 0; i < res.size(); i++){
        for (int j = 0; j < res[0].size(); j++){
            float_res[i][j] = ((float)int_res[i][j])/((float)pow(2,length/2));
        }
    }

    cout << "### Output Data ###" << endl;
    for (int row = 0; row < res.size(); row++) {
        for (int col = 0; col < res[0].size(); col++) {
            cout << float_res[row][col] << "  ";
        }
        cout << "\n";
    }
    cout << "\n";

    delete_gate_bootstrapping_secret_keyset(key);
}


vector<vector<float>> DomEigenvector(vector<vector<float>> matrix, const int number_of_iteration) {
    vector<vector<float>> result(matrix.size(), vector<float>(1));
    vector<vector<float>> result_temp(matrix.size(), vector<float>(1));
    vector<vector<float>> res(matrix.size(), vector<float>(1));
    vector<vector<float>> plainres(matrix.size(), vector<float>(1));
    for (int i = 0; i < matrix.size(); i++){
        result[i][0] = 0;
    }
    for (int i = 0; i < matrix.size(); i++){
        for (int j = 0 ; j < matrix.size(); j++){
            result[i][0] += matrix[i][j];
        }
    }
    for (int i = 0; i < matrix.size(); i++){
        result_temp[i][0] = abs(result[i][0]);
    }
    for (int i = 0; i < matrix.size()-1; i++){
        if(result_temp[i][0] > result_temp[i+1][0]){
            result_temp[i+1][0] = result_temp[i][0];
        }
    }
    float max_value = result_temp[matrix.size()-1][0];
    float norm_constant = (float)1/max_value;
    for (int i = 0; i < matrix.size(); i++){
            res[i][0] = result[i][0]*norm_constant;
    }

    cout << "after 1st iteration\n";
    for(int i = 0; i < matrix.size(); i++){
        cout << res[i][0] << "\n";
    }
    cout << "\n";

    for (int h = 0; h < number_of_iteration - 1; h++){
        for (int i = 0; i < matrix.size(); i++){
            result[i][0] = 0;
        }
        for(int i = 0; i < matrix.size(); i++){
            for(int j = 0; j < 1; j++){
                for(int k = 0; k < matrix.size(); k++){
                    result[i][j] += matrix[i][k] * res[k][j];
                }
            }
        }
        for (int i = 0; i < matrix.size(); i++){
            result_temp[i][0] = abs(result[i][0]);
        }
        for (int i = 0; i < matrix.size()-1; i++){
            if(result_temp[i][0] > result_temp[i+1][0]){
                result_temp[i+1][0] = result_temp[i][0];
            }
        }
        float max_value = result_temp[matrix.size()-1][0];
        float norm_constant = (float)1/max_value;
        for (int i = 0; i < matrix.size(); i++){
                res[i][0] = result[i][0]*norm_constant;
        }
    
        cout << "after %dth iteration\n", h+2;
        for(int i = 0; i < matrix.size(); i++){
            cout << res[i][0] << "\n";
        }
        cout << "\n";
    }
    for(int i = 0; i < matrix.size(); i++){
        plainres[i][0] = res[i][0];
    }

    return plainres;
}


// res = nx1 matrix, matrix = nxn matrix
// QR method
vector<vector<float>> DomEigenvector_QR(vector<vector<float>> matrix, const int number_of_iteration) {
    vector<vector<float>> epsil(matrix.size(), vector<float>(1));
    vector<vector<float>> A(matrix.size(), vector<float>(1));
    vector<vector<float>> B(matrix.size(), vector<float>(1));
    vector<vector<float>> transpose_B(1, vector<float>(matrix.size()));
    vector<vector<float>> mult_RB(matrix.size(), vector<float>(1));
    vector<vector<float>> R(matrix[0].size(), vector<float>(matrix[0].size()));
    vector<vector<float>> Q(matrix.size(), vector<float>(matrix[0].size()));
    vector<vector<float>> Matrix(matrix.size(), vector<float>(matrix[0].size()));
    vector<vector<float>> EigenvectorMatrix(matrix.size(), vector<float>(matrix[0].size()));
    vector<vector<float>> r(1, vector<float>(1));
    vector<vector<float>> plainres(matrix.size(), vector<float>(matrix[0].size()));
    
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            Matrix[i][j] = matrix[i][j];
            EigenvectorMatrix[i][j] = 0;
        }
        EigenvectorMatrix[i][i] = 1;
    }

    for(int n = 0; n < number_of_iteration; n++){
        /*QR decomposition Matrix=QR*/
        for(int i = 0; i < matrix.size(); i++){
            epsil[i][0] = Matrix[i][0];
        }
        for(int i = 0; i < matrix[0].size(); i++){
            for(int j = 0; j < matrix[0].size(); j++){
                R[i][j] = 0;
            }   
        }
        float norm_epsil = 0;

        for(int i = 0; i < matrix.size(); i++){
            norm_epsil += (epsil[i][0])*(epsil[i][0]);
        }
        
        R[0][0] = norm_epsil;

        for(int i = 0; i < matrix.size(); i++){
            Q[i][0] = (epsil[i][0])*(1/R[0][0]);
        }

        for(int j = 1; j < matrix[0].size(); j++){
            for(int i = 0; i < matrix.size(); i++){
                A[i][0] = Matrix[i][j];
                epsil[i][0] = Matrix[i][j];
            }

            for(int h = 0; h < j; h++){
                for(int i = 0; i < matrix.size(); i++){
                    B[i][0] = Q[i][h];
                }
                transpose_B = TransposeMatrix(B);
                r = MultiplyMatrix(transpose_B, A);
                R[h][j] = r[0][0];
                mult_RB = MultiplyMatrix(B, r);
                epsil = SubtractMatrix(epsil, mult_RB);
            }
        
            norm_epsil = 0;
            for(int i = 0; i < matrix.size(); i++){
                norm_epsil += (epsil[i][0])*(epsil[i][0]);
            }
            R[j][j] = norm_epsil;

            for(int i = 0; i < matrix.size(); i++){
                Q[i][j] = (epsil[i][0])*(1/norm_epsil);
            }
        }

        /*Matrix=RQ, V=VQ update*/
        Matrix = MultiplyMatrix(R, Q);
        EigenvectorMatrix = MultiplyMatrix(EigenvectorMatrix, Q);
    }

    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            plainres[i][j] = EigenvectorMatrix[i][j];
        }
    }

    return plainres;
}





void newprint_Matrix(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    for (int i = 0; i < matrix.size(); i++) { 
       	for(int j = 0; j < matrix[0].size(); j++){
            for(int k = 0; k < length; k++){
                bootsCOPY(&res[i][j][k], &matrix[i][j][k], bk);
            }
        }
	}
}


void newAddMatrix(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, vector<vector<LweSample*>> bMatrix, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    for (int row = 0; row < aMatrix.size(); row++) {
        for (int col = 0; col < aMatrix[0].size(); col++) {
            newADD(res[row][col], aMatrix[row][col], bMatrix[row][col], length, bk);
        }
    }
}


void newmultRealMatrix(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, vector<vector<LweSample*>> bMatrix, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    vector<LweSample*> temp(aMatrix[0].size());
    for(int i = 0; i < aMatrix[0].size(); i++){
        temp[i] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    }
    vector<vector<LweSample*>> A(aMatrix.size(), vector<LweSample*> (aMatrix[0].size()));
    for(int i = 0; i < aMatrix.size(); i++){
        for(int j = 0; j < aMatrix[0].size(); j++){
            A[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
    }
    vector<vector<LweSample*>> B(bMatrix.size(), vector<LweSample*> (bMatrix[0].size()));
    for(int i = 0; i < bMatrix.size(); i++){
        for(int j = 0; j < bMatrix[0].size(); j++){
            B[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
    }

    for(int i = 0; i < aMatrix.size(); i++){
        for(int j = 0; j < aMatrix[0].size(); j++){
            for(int k = 0; k < length; k++){
                bootsCOPY(&A[i][j][k], &aMatrix[i][j][k], bk);
            }
        }
    }
    for(int i = 0; i < bMatrix.size(); i++){
        for(int j = 0; j < bMatrix[0].size(); j++){
            for(int k = 0; k < length; k++){
                bootsCOPY(&B[i][j][k], &bMatrix[i][j][k], bk);
            }
        }
    }

    for (int i = 0; i < aMatrix.size(); i++) {
        for (int j = 0; j < bMatrix[0].size(); j++) {
            for (int k = 0; k < length; k++){
                bootsCONSTANT(&res[i][j][k], 0, bk);
            }
        }
    }

    // perform matrix multiplication
    for (int i = 0; i < aMatrix.size(); i++) {
        for (int j = 0; j < bMatrix[0].size(); j++) {
            #pragma omp parallel for
            for (int k = 0; k < aMatrix[0].size(); k++) {
                newMultiReal(temp[k], A[i][k], B[k][j], length, bk);               
            }
     
            for (int k = 0; k < aMatrix[0].size(); k++) {
                newADD(res[i][j], res[i][j], temp[k], length, bk);
            }
        }
    }
    
    for(int i = 0; i < aMatrix[0].size(); i++){
        delete_gate_bootstrapping_ciphertext_array(length, temp[i]);
    }
    for(int i = 0; i < aMatrix.size(); i++){
        for(int j = 0; j < aMatrix[0].size(); j++){
            delete_gate_bootstrapping_ciphertext_array(length, A[i][j]);
        }
    }
    for(int i = 0; i < bMatrix.size(); i++){
        for(int j = 0; j < bMatrix[0].size(); j++){
            delete_gate_bootstrapping_ciphertext_array(length, B[i][j]);
        }
    }
}


// void newmultRealMatrix2(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, vector<vector<LweSample*>> bMatrix, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
// 	LweSample* result[aMatrix.size()][bMatrix[0].size()];
// 	for(int i = 0; i < aMatrix.size(); i++){
//         for(int j = 0; j < bMatrix[0].size(); j++){
// 		    result[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
//         }
// 	}        
//     LweSample* temp = new_gate_bootstrapping_ciphertext_array(length, bk->params);

//     for (int i = 0; i < aMatrix.size(); i++) {
//         for (int j = 0; j < bMatrix[0].size(); j++) {
//             bootsCONSTANT(result[i][j], 0, bk);
//             bootsCONSTANT(res[i][j], 0, bk);
//         }
//     }

//     // perform matrix multiplication
//     for (int i = 0; i < aMatrix.size(); i++) {
//         for (int j = 0; j < bMatrix[0].size(); j++) {
//         //     for (int k = 0; k < M; k++) {
//         //         newMultiReal2(temp, a[i][k], b[k][j], length, bk);
//         //         newADD(result[i][j], result[i][j], temp, length, bk);
//         //     } 
//         newMultiReal2(result[i][j], aMatrix[i][j], bMatrix[i][j], length, bk);
//         }
//     }
    
//     for (int i = 0; i < aMatrix.size(); i++) { 
//        	for(int j = 0; j < bMatrix[0].size(); j++){
//             for(int k = 0; k < length; k++){
// 		        bootsCOPY(&res[i][j][k], &result[i][j][k], bk);
//             }
//         }
// 	}

//     delete_gate_bootstrapping_ciphertext_array(length, temp);
// 	for(int i = 0; i < aMatrix.size(); i++){
//         for(int j = 0; j < bMatrix[0].size(); j++){
// 		    delete_gate_bootstrapping_ciphertext_array(length, result[i][j]);
//         }
// 	}
// }


void newSubMatrix(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, vector<vector<LweSample*>> bMatrix, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    for (int row = 0; row < aMatrix.size(); row++) {
        for (int col = 0; col < aMatrix[0].size(); col++) {
            newSUB(res[row][col], aMatrix[row][col], bMatrix[row][col], length, bk);
        }
    }
}


void newTraceMatrix(LweSample* res, vector<vector<LweSample*>> matrix, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    // sum of diag element
    for (int i = 0; i < length; i++) {
        bootsCONSTANT(&res[i], 0, bk);
    }

    // Add diag element
    for (int i = 0; i < matrix.size(); i++) {
        newADD(res, res, matrix[i][i], length, bk);
    }
}


void newTransposeMatrix(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    for (int row = 0; row < matrix.size(); row++) {
        for (int col = 0; col < matrix[0].size(); col++) {
            for (int k = 0; k < length; k++){
                bootsCOPY(&res[col][row][k], &matrix[row][col][k], bk);
            }
        }
    }
}


void newP2C(LweSample* res, const float num, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    int32_t plain = (int32_t)((num)*pow(2,length/2));
    // for(int i = 0; i < (length/2); i++){
    //     bootsCONSTANT(&res[i], 0, bk);
    // }
    // for(int i = 0; i < (length/2); i++){
    //     bootsCONSTANT(&res[i+(length/2)], (plain>>i)&1, bk);
    // }
    for(int i = 0; i < length; i++){
        bootsCONSTANT(&res[i], (plain>>i)&1, bk);
    }
}


void newScalarMultMatrix(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, const float num, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    LweSample* EN = new_gate_bootstrapping_ciphertext_array(length, bk->params);

    newP2C(EN, num, length, bk);

    for (int i = 0; i < matrix.size(); i++){
        for (int j = 0; j < matrix[0].size(); j++){
            newMultiReal(res[i][j], matrix[i][j], EN, length, bk);
        }
    }

    delete_gate_bootstrapping_ciphertext_array(length, EN);
}


void newScalarMultMatrix2(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[0].size(); j++) {
            newMultiReal(res[i][j], matrix[i][j], b, length, bk);
        }
    }
}


void newInverseMatrix(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, const int number_of_iteration, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    vector<vector<LweSample*>> transposed_a(matrix[0].size(), vector<LweSample*>(matrix.size()));
    vector<vector<LweSample*>> result(matrix[0].size(), vector<LweSample*>(matrix.size()));
    for(int i = 0; i < matrix[0].size(); i++){
        for(int j = 0; j < matrix.size(); j++){
            transposed_a[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
            result[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
    }

    vector<vector<LweSample*>> matrix_multipled_a(matrix.size(), vector<LweSample*>(matrix.size()));
    vector<vector<LweSample*>> identity_matrix(matrix.size(), vector<LweSample*>(matrix.size()));
    vector<vector<LweSample*>> ab(matrix.size(), vector<LweSample*>(matrix.size()));
    vector<vector<LweSample*>> temp(matrix.size(), vector<LweSample*>(matrix.size()));
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix.size(); j++){
            matrix_multipled_a[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
            identity_matrix[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
            ab[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
            temp[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
    }

    LweSample* E1 = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* trace = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* alpha = new_gate_bootstrapping_ciphertext_array(length, bk->params);


    // Calculate aa^T
    newTransposeMatrix(transposed_a, matrix, length, bk);
    newmultRealMatrix(matrix_multipled_a, matrix, transposed_a, length, bk); 


    // Find Trace = Diag element sum
    newTraceMatrix(trace, matrix_multipled_a, length, bk);

    // Let 1/Trace = alpha. Calculate alpha * transposed_a
    float num = 1;
    newP2C(E1, num, length, bk); 
    newRealDiv(alpha, E1, trace, length, bk); 
    newScalarMultMatrix2(result, transposed_a, alpha, length, bk);  

    // Set the identity matrix - Instead multipy 2 to identity_matrix, putting 1 to 17-bit
    for (int i = 0; i < matrix.size(); i++){
        for (int j = 0; j < matrix.size(); j++){
            for (int k = 0; k < length; k++){
                bootsCONSTANT(&identity_matrix[i][j][k], 0, bk);
            }
        }
    }
    /*for (int i = 0; i < matrix.size(); i++){
        bootsCONSTANT(&identity_matrix[i][i][16], 1, bk); 
    }*/

    for (int i = 0; i < number_of_iteration; i++){
        newmultRealMatrix(ab, matrix, result, length, bk);  
        newSubMatrix(temp, identity_matrix, ab, length, bk);
        newmultRealMatrix(res, result, temp, length, bk); 
    }


    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix.size(); j++){
            delete_gate_bootstrapping_ciphertext_array(length, matrix_multipled_a[i][j]);
            delete_gate_bootstrapping_ciphertext_array(length, identity_matrix[i][j]);
            delete_gate_bootstrapping_ciphertext_array(length, ab[i][j]);
            delete_gate_bootstrapping_ciphertext_array(length, temp[i][j]);
        }
    }
    for(int i = 0; i < matrix[0].size(); i++){
        for(int j = 0; j < matrix.size(); j++){
            delete_gate_bootstrapping_ciphertext_array(length, transposed_a[i][j]);
            delete_gate_bootstrapping_ciphertext_array(length, result[i][j]);
        }
    }
    delete_gate_bootstrapping_ciphertext_array(length, alpha);
    delete_gate_bootstrapping_ciphertext_array(length, trace);
    delete_gate_bootstrapping_ciphertext_array(length, E1);
}

//res = nx1 matrix, matrix = nxn matrix
void newDomEigenvector(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, const int number_of_iteration, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    vector<vector<LweSample*>> result(matrix.size(), vector<LweSample*>(1));
    vector<vector<LweSample*>> result_temp(matrix.size(), vector<LweSample*>(1));
    for(int i = 0; i < matrix.size(); i++){
        result[i][0] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        result_temp[i][0] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	}
    LweSample* max_value = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* E1 = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* norm_constant = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    
    // Set an initial nonzero N*1 vector(matrix) - we make its every element 1
    // Instead of matrix multiple, we use adding - since we used the N*1 matrix with every element 1
    for (int i = 0; i < matrix.size(); i++){
        for (int j = 0 ; j < length; j++){
            bootsCONSTANT(&result[i][0][j], 0, bk);
        }
    }

    for (int i = 0; i < matrix.size(); i++){
        for (int j = 0; j < matrix.size(); j++){
            newADD(result[i][0], result[i][0], matrix[i][j], length, bk);
        }
    }

    // cout << "domeigenvector _ 1st result" << endl;
    // decrypting_code(result, length);

    for (int j = 0; j < matrix.size(); j++){
        newABS(result_temp[j][0], result[j][0], length, bk);
    }
    
    // cout << "domeigenvector _ 1st result_temp" << endl;
    // decrypting_code(result_temp, length);

    newMaxValue(max_value, result_temp, length, bk);    
    newP2C(E1, 1, length, bk);    
   
    newRealDiv(norm_constant, E1, max_value, length, bk);
    
    for(int i= 0; i < length; i++){
        bootsCOPY(&res[0][0][i], &norm_constant[i], bk);
    }
    // cout << "domeigenvector _ 1st norm_constant" << endl;
    // decrypting_code(res, length);

    newScalarMultMatrix2(res, result, norm_constant, length, bk);
     
    cout << "domeigenvector _ after 1st iteration\n";
    decrypting_code(res, length);

    for (int i = 0; i < number_of_iteration - 1; i++){
        newmultRealMatrix(result, matrix, res, length, bk);
        // Normalize the vectors - Divide the vectors by the maximum value of absolute number of the elements
      
        for (int j = 0; j < matrix.size(); j++){
            newABS(result_temp[j][0], result[j][0], length, bk);
      
        }

        // cout << "domeigenvector _ %dth result_temp", i+2 << endl;
        // decrypting_code(result_temp, length);

        newMaxValue(max_value, result_temp, length, bk);  
        newRealDiv(norm_constant, E1, max_value, length, bk);
        newScalarMultMatrix2(res, result, norm_constant, length, bk);

        cout << "domeigenvector _ after %dth iteration\n", i+2 << endl;
        decrypting_code(res, length);
    }


    for(int i = 0; i < matrix.size(); i++){
	    delete_gate_bootstrapping_ciphertext_array(length, result[i][0]);
        delete_gate_bootstrapping_ciphertext_array(length, result_temp[i][0]);
	}
    delete_gate_bootstrapping_ciphertext_array(length, max_value);
    delete_gate_bootstrapping_ciphertext_array(length, E1);
    delete_gate_bootstrapping_ciphertext_array(length, norm_constant);
}


// Instead of using Rayleigh quotient, we use matrix multiple of first element of eigenvector.
// x[0] / x_prime[0] = lambda. ( x_prime is the vector with number_of_iteration - 1 / [0] means the first element of the matrix.)
void newDomEigenvalue(LweSample* res, vector<vector<LweSample*>> matrix, const int number_of_iteration, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    vector<vector<LweSample*>> eigenvector_prime(matrix.size(), vector<LweSample*> (1));
    vector<vector<LweSample*>> eigenvector(matrix.size(), vector<LweSample*> (1));
    for(int i = 0; i < matrix.size(); i++){
        eigenvector[i][0] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        eigenvector_prime[i][0] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    }
    LweSample* E1 = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* denominator = new_gate_bootstrapping_ciphertext_array(length, bk->params);

    newP2C(E1, 1, length, bk);

    // Find the x_prime with number_of_iteration-1.
    newDomEigenvector(eigenvector_prime, matrix, number_of_iteration-1, length, bk);

    // With no normalization, just find the x[0].
    newmultRealMatrix(eigenvector, matrix, eigenvector_prime, length, bk);

    // Divide x[0] by x_prime[0].
    newRealDiv(denominator, E1, eigenvector_prime[0][0], length, bk);
    newMultiReal(res, eigenvector[0][0], denominator, length, bk);

    for(int i = 0; i < matrix.size(); i++){
        delete_gate_bootstrapping_ciphertext_array(length, eigenvector[i][0]);
        delete_gate_bootstrapping_ciphertext_array(length, eigenvector_prime[i][0]);
    }
    delete_gate_bootstrapping_ciphertext_array(length, E1);
    delete_gate_bootstrapping_ciphertext_array(length, denominator);
}



// res = nx1 matrix, matrix = nxn matrix
// QR method
void newDomEigenvector_QR(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, const int number_of_iteration, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    vector<vector<LweSample*>> epsil(matrix.size(), vector<LweSample*>(1));
    vector<vector<LweSample*>> A(matrix.size(), vector<LweSample*>(1));
    vector<vector<LweSample*>> B(matrix.size(), vector<LweSample*>(1));
    vector<vector<LweSample*>> transpose_B(1, vector<LweSample*>(matrix.size()));
    vector<vector<LweSample*>> mult_RB(matrix.size(), vector<LweSample*>(1));
    for(int i = 0; i < matrix.size(); i++){
        epsil[i][0] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        A[i][0] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        B[i][0] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        transpose_B[0][i] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        mult_RB[i][0] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    }
    vector<vector<LweSample*>> R(matrix[0].size(), vector<LweSample*>(matrix[0].size()));
    for(int i = 0; i < matrix[0].size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            R[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
    }
    vector<vector<LweSample*>> Q(matrix.size(), vector<LweSample*>(matrix[0].size()));
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            Q[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
    }
    LweSample* P2C_1 = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* P2C_1_norm = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* square = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* norm_epsil = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    vector<vector<LweSample*>> Matrix(matrix.size(), vector<LweSample*>(matrix[0].size()));
    vector<vector<LweSample*>> EigenvectorMatrix(matrix.size(), vector<LweSample*>(matrix[0].size()));
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            Matrix[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
            EigenvectorMatrix[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
    }
    vector<vector<LweSample*>> r(1, vector<LweSample*>(1));
    r[0][0] = new_gate_bootstrapping_ciphertext_array(length, bk->params);


    newP2C(P2C_1, 1, length, bk);
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            for(int k = 0; k < length; k++){
                bootsCOPY(&Matrix[i][j][k], &matrix[i][j][k], bk);
                bootsCONSTANT(&EigenvectorMatrix[i][j][k], 0, bk);
            }
        }
        bootsCONSTANT(&EigenvectorMatrix[i][i][length/2], 1, bk);
    }

    for(int n = 0; n < number_of_iteration; n++){
        /*QR decomposition Matrix=QR*/
        for(int i = 0; i < matrix.size(); i++){
            for(int k = 0; k < length; k++){
                bootsCOPY(&epsil[i][0][k], &Matrix[i][0][k], bk);
            }
        }
        for(int i = 0; i < matrix[0].size(); i++){
            for(int j = 0; j < matrix[0].size(); j++){
                for(int k = 0; k < length; k++){
                    bootsCONSTANT(&R[i][j][k], 0, bk);
                }
            }
        }
        for(int k = 0; k < length; k++){
            bootsCONSTANT(&norm_epsil[k], 0, bk);
        }

        for(int i = 0; i < matrix.size(); i++){
            newMultiReal(square, epsil[i][0], epsil[i][0], length, bk);
            newADD(norm_epsil, norm_epsil, square, length, bk);             //||epsil_0||
        }
        for(int k = 0; k < length; k++){
            bootsCOPY(&R[0][0][k], &norm_epsil[k], bk);
        }

        newRealDiv(P2C_1_norm, P2C_1, norm_epsil, length, bk);
        for(int i = 0; i < matrix.size(); i++){
            newMultiReal(Q[i][0], epsil[i][0], P2C_1_norm, length, bk);     //Q_0
        }

        for(int j = 1; j < matrix[0].size(); j++){
            for(int i = 0; i < matrix.size(); i++){
                for(int k = 0; k < length; k++){
                    bootsCOPY(&A[i][0][k], &Matrix[i][j][k], bk);           //A_1
                    bootsCOPY(&epsil[i][0][k], &Matrix[i][j][k], bk);
                }
            }

            for(int h = 0; h < j; h++){
                for(int i = 0; i < matrix.size(); i++){
                    for(int k = 0; k < length; k++){
                        bootsCOPY(&B[i][0][k], &Q[i][h][k], bk);            //B = Q_0
                    }
                }
                newTransposeMatrix(transpose_B, B, length, bk);
                newmultRealMatrix(r, transpose_B, A, length, bk);
                for(int k = 0; k < length; k++){
                    bootsCOPY(&R[h][j][k], &r[0][0][k], bk);
                }
                newScalarMultMatrix2(mult_RB, B, R[h][j], length, bk);
                newSubMatrix(epsil, epsil, mult_RB, length, bk);            //epsil_1
            }

            for(int k = 0; k < length; k++){
                bootsCONSTANT(&norm_epsil[k], 0, bk);
            }
            for(int i = 0; i < matrix.size(); i++){
                newMultiReal(square, epsil[i][0], epsil[i][0], length, bk);
                newADD(norm_epsil, norm_epsil, square, length, bk);         //||epsil_1||
            }
            for(int k = 0; k < length; k++){
                bootsCOPY(&R[j][j][k], &norm_epsil[k], bk);
            }

            newRealDiv(P2C_1_norm, P2C_1, norm_epsil, length, bk);
            for(int i = 0; i < matrix.size(); i++){
                newMultiReal(Q[i][j], epsil[i][0], P2C_1_norm, length, bk); //Q_1
            }
        }

        /*Matrix=RQ, V=VQ update*/
        newmultRealMatrix(Matrix, R, Q, length, bk);
        newmultRealMatrix(EigenvectorMatrix, EigenvectorMatrix, Q, length, bk);
    }

    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            for(int k = 0; k < length; k++){
                bootsCOPY(&res[i][j][k], &EigenvectorMatrix[i][j][k], bk);
            }
        }
    }


    for(int i = 0; i < matrix.size(); i++){
        delete_gate_bootstrapping_ciphertext_array(length, epsil[i][0]);
        delete_gate_bootstrapping_ciphertext_array(length, A[i][0]);
        delete_gate_bootstrapping_ciphertext_array(length, B[i][0]);
        for(int j = 0; j < matrix[0].size(); j++){
            delete_gate_bootstrapping_ciphertext_array(length, Q[i][j]);
        }
        delete_gate_bootstrapping_ciphertext_array(length, transpose_B[0][i]);
        delete_gate_bootstrapping_ciphertext_array(length, mult_RB[i][0]);
    }
    for(int i = 0; i < matrix[0].size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            delete_gate_bootstrapping_ciphertext_array(length, R[i][j]);
        }
    }
    delete_gate_bootstrapping_ciphertext_array(length, P2C_1);
    delete_gate_bootstrapping_ciphertext_array(length, P2C_1_norm);
    delete_gate_bootstrapping_ciphertext_array(length, square);
    delete_gate_bootstrapping_ciphertext_array(length, norm_epsil);
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            delete_gate_bootstrapping_ciphertext_array(length, Matrix[i][j]);
            delete_gate_bootstrapping_ciphertext_array(length, EigenvectorMatrix[i][j]);
        }
    }
    delete_gate_bootstrapping_ciphertext_array(length, r[0][0]);
}




// function getting r eigenvectors
// qr method
// void {
//     qr method
//     -> eigenvalue
//     -> pick eigenvectors according to eigenvalue size
// }


// function getting r eigenvectors
// power method
// void {
//     power iteration method eigenvector
//     -> r times
//     -> make it into matrix. done.
// }




void practice(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, vector<vector<LweSample*>> bMatrix, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    vector<vector<LweSample*>> temp(aMatrix.size(), vector<LweSample*>(aMatrix[0].size()));        
	for(int i = 0; i < aMatrix.size(); i++){
        for(int j = 0; j < aMatrix[0].size(); j++){
            temp[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
	}
    newAddMatrix(temp, aMatrix, bMatrix, length, bk);
    decrypting_code(temp, length);
    newAddMatrix(temp, temp, temp, length, bk);
    decrypting_code(temp, length);
    for(int i = 0; i < temp.size(); i++){
        for(int j = 0; j < temp[0].size(); j++){
            for(int k = 0; k < length; k++){
                bootsCOPY(&res[i][j][k], &temp[i][j][k], bk);
            }
        }
    }
    decrypting_code(res, length);
    newmultRealMatrix(res, res, aMatrix, length, bk);
    decrypting_code(res, length);
    
    newMultiReal(res[0][0], res[0][0], temp[0][0], length, bk);
}
    