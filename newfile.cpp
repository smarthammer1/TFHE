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
    for (int h = 0; h < number_of_iteration - 1; h++){
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
    }
    for(int i = 0; i < matrix.size(); i++){
        plainres[i][0] = res[i][0];
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
    // LweSample* temp = new_gate_bootstrapping_ciphertext_array(length, bk->params);

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
                newMultiReal(temp[k], aMatrix[i][k], bMatrix[k][j], length, bk);               
            } 
            for (int k = 0; k < aMatrix[0].size(); k++) {
                newADD(res[i][j], res[i][j], temp[k], length, bk);
            }
        }
    }
    
    for(int i = 0; i < aMatrix[0].size(); i++){
        delete_gate_bootstrapping_ciphertext_array(length, temp[i]);
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
    for (int j = 0; j < matrix.size(); j++){
        newABS(result_temp[j][0], result[j][0], length, bk);
    }
    newMaxValue(max_value, result_temp, length, bk);    
    newP2C(E1, 1, length, bk);    
    newRealDiv(norm_constant, E1, max_value, length, bk);
    newScalarMultMatrix2(res, result, norm_constant, length, bk);
    for (int i = 0; i < number_of_iteration - 1; i++){
        newmultRealMatrix(result, matrix, res, length, bk);
        // Normalize the vectors - Divide the vectors by the maximum value of absolute number of the elements
        for (int j = 0; j < matrix.size(); j++){
            newABS(result_temp[j][0], result[j][0], length, bk);
        }
        newMaxValue(max_value, result_temp, length, bk);  
        newRealDiv(norm_constant, E1, max_value, length, bk);
        newScalarMultMatrix2(res, result, norm_constant, length, bk);
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


