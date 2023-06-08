#include <cmath>
using namespace std;

vector<float> PCAMatrix_step1(vector<vector<float>> matrix) {
	vector<float> plainres(matrix[0].size());
	vector<float> centering(matrix[0].size());
	for(int i = 0; i < matrix[0].size(); i++){
		centering[i] = 0;
	}
    float E1_EN = (float)1/matrix.size();
    for(int j = 0; j < matrix[0].size(); j++){
        for(int i = 0; i < matrix.size(); i++){
            centering[j] += matrix[i][j];
        }
    }
    for(int j = 0; j < matrix[0].size(); j++){
        plainres[j] = (centering[j])*(E1_EN);
    }
    return plainres;
}

vector<vector<float>> PCAMatrix_step2(vector<vector<float>> matrix, vector<float> centering) {
    vector<vector<float>> plainres(matrix[0].size(), vector<float> (matrix[0].size()));
    vector<vector<float>> mean_centered_matrix(matrix.size(), vector<float> (matrix[0].size()));
    vector<vector<float>> transpose_mean_centered_matrix(matrix[0].size(), vector<float> (matrix.size()));
    vector<vector<float>> S(matrix[0].size(), vector<float> (matrix[0].size()));
    for(int i = 0; i < matrix[0].size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            S[i][j] = 0;
        }
    }

    float E1_EN = (float)1/matrix.size();
    
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            mean_centered_matrix[i][j] = matrix[i][j] - centering[j];
        }
    }
    for(int i = 0; i < matrix[0].size(); i++){
        for(int j = 0; j < matrix.size(); j++){
            transpose_mean_centered_matrix[i][j] = mean_centered_matrix[j][i];
        }
    }
    for(int i = 0; i < matrix[0].size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            for(int k = 0; k < matrix.size(); k++){
                S[i][j] += transpose_mean_centered_matrix[i][k] * mean_centered_matrix[k][j];
            }
            plainres[i][j] = E1_EN * S[i][j];
        }
    }
    return plainres;
}


vector<vector<float>> PCAMatrix_step4(vector<vector<float>> eigenvalue_S, vector<vector<float>> eigenvector_S, const int reduced_basis) {
    vector<vector<float>> plainres(eigenvector_S.size(), vector<float> (reduced_basis));    
    vector<float> temp(2);
    for(int i = 0; i < 2; i++){
        temp[i] = 0;
    }
    float tmp, a, b;

    for (int h = 0; h < eigenvalue_S.size() - 1; h++) {
        for (int j = 0; j < eigenvalue_S.size() - 1 - h; j++) {
            if(eigenvalue_S[j][0] < eigenvalue_S[j+1][0]){
                a = eigenvalue_S[j][0];
                eigenvalue_S[j][0] = eigenvalue_S[j+1][0];
                eigenvalue_S[j+1][0] = a;
                for(int i = 0; i < eigenvector_S.size(); i++){
                    b = eigenvector_S[i][j];
                    eigenvector_S[i][j] = eigenvector_S[i][j+1];
                    eigenvector_S[i][j+1] = b;
                }
            }
        }
    }
    for(int i = 0; i < eigenvector_S.size(); i++){
        for(int j = 0; j < reduced_basis; j++){
            plainres[i][j] = eigenvector_S[i][j];
        }
    }
    return plainres;
}


vector<vector<float>> PCAMatrix_step5(vector<vector<float>> aMatrix, vector<vector<float>> reduced_Matrix) {
    vector<vector<float>> plainres(reduced_Matrix[0].size(), vector<float>(aMatrix.size()));
    vector<vector<float>> transpose_reducedMatrix(reduced_Matrix[0].size(), vector<float>(reduced_Matrix.size()));
    vector<vector<float>> transpose_aMatrix(aMatrix[0].size(), vector<float>(aMatrix.size()));
    vector<vector<float>> res(reduced_Matrix[0].size(), vector<float>(aMatrix.size()));

    // Project the data onto the new subspace by computing Z' = ZW   // compute reduced dimensionality data
    transpose_reducedMatrix = TransposeMatrix(reduced_Matrix);
    transpose_aMatrix = TransposeMatrix(aMatrix);
    res = MultiplyMatrix(transpose_reducedMatrix, transpose_aMatrix);

    for (int row = 0; row < reduced_Matrix[0].size(); row++) {
        for (int col = 0; col < aMatrix.size(); col++) {
            plainres[row][col] =  res[row][col];
        }
    }
    return plainres;
}


vector<vector<float>> PCAMatrix_step3_4_reducedbasis(vector<vector<float>> matrix, const int number_of_iteration, const int dimension) {
    vector<vector<float>> plainres(matrix.size(), vector<float>(dimension));  
    vector<vector<float>> result(matrix.size(), vector<float>(1));
    vector<vector<float>> eigenvector(matrix.size(), vector<float> (1));
    vector<vector<float>> eigenvector2(matrix.size(), vector<float> (1));
    float num;
    vector<vector<float>> tmp(matrix.size(), vector<float>(matrix.size()));
    vector<vector<float>> reduced_basisMatrix(matrix.size(), vector<float>(dimension));        
    float sqsum;
    float square;
    float sqroot;
    vector<vector<float>> normalized(matrix.size(), vector<float>(1));
    vector<vector<float>> transpose_normalized(1, vector<float>(matrix.size()));

    for(int h = 0; h < dimension - 1; h++){    
        result = DomEigenvector(matrix, number_of_iteration - 1);
        eigenvector = DomEigenvector(matrix, number_of_iteration);
        eigenvector2 = MultiplyMatrix(matrix, result);                    /////////this is the dominant eigenvector of S. doing one more iteration on result.

        float num = eigenvector2[0][0] / result[0][0];       ///eigenvalue end///// eigenvalue = num

        float sqsum = 0;
        for(int i = 0; i < matrix.size(); i++){
            sqsum += eigenvector[i][0]*eigenvector[i][0];
        }
        float sqroot = sqrt(sqsum);                         // sqareroot of sqsum

        for(int i = 0; i < matrix.size(); i++){
            normalized[i][0] = (1/sqroot)*eigenvector[i][0];         
        }
        transpose_normalized = TransposeMatrix(normalized);

        tmp = MultiplyMatrix(normalized, transpose_normalized);
        for(int i = 0; i < matrix.size(); i++){
            for(int j = 0; j < matrix.size(); j++){
                tmp[i][j] = num*tmp[i][j];
            }
        }

        matrix = SubtractMatrix(matrix, tmp);                            ///////deflated matrix

        for(int i = 0; i < matrix.size(); i++){
            reduced_basisMatrix[i][h] = eigenvector[i][0];
        }
    }

    result = DomEigenvector(matrix, number_of_iteration);
    for(int i = 0; i < matrix.size(); i++){
        reduced_basisMatrix[i][dimension - 1] = result[i][0];
    }


    for (int row = 0; row < matrix.size(); row++) {
        for (int col = 0; col < dimension; col++) {
            plainres[row][col] =  reduced_basisMatrix[row][col];
        }
    }
    return plainres;
}




// // change number of dimensionality -> then change everything else related to it.
vector<vector<float>> PCAMatrix(vector<vector<float>> matrix, const int number_of_iteration, const int dimension){
    vector<vector<float>> plainres(matrix.size(), vector<float>(dimension));    
    vector<float> centering(matrix[0].size());
    vector<vector<float>> mean_centered_matrix(matrix.size(), vector<float> (matrix[0].size()));
    vector<vector<float>> transpose_mean_centered_matrix(matrix[0].size(), vector<float> (matrix.size()));
    vector<vector<float>> S(matrix[0].size(), vector<float> (matrix[0].size()));
    vector<vector<float>> normalized(matrix[0].size(), vector<float>(1));
    vector<vector<float>> transpose_normalized(1, vector<float>(matrix[0].size()));
    vector<vector<float>> tmp(matrix[0].size(), vector<float>(matrix[0].size()));
    vector<vector<float>> result(matrix[0].size(), vector<float>(1));
    vector<vector<float>> reduced_basisMatrix(matrix[0].size(), vector<float>(dimension));  
    vector<vector<float>> eigenvector(matrix[0].size(), vector<float>(1));
    vector<vector<float>> eigenvector2(matrix[0].size(), vector<float>(1));
    vector<vector<float>> result_temp(matrix[0].size(), vector<float>(1));

    // step1
	for(int i = 0; i < matrix[0].size(); i++){
		centering[i] = 0;
	}
    for(int j = 0; j < matrix[0].size(); j++){
        for(int i = 0; i < matrix.size(); i++){
            centering[j] += matrix[i][j];
        }
    }
    float E1_EN = (float)1/matrix.size();
    for(int j = 0; j < matrix[0].size(); j++){
        centering[j] = (centering[j])*(E1_EN);                                             
    }

    // step2
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            mean_centered_matrix[i][j] = matrix[i][j] - centering[j];
        }
    }
    transpose_mean_centered_matrix = TransposeMatrix(mean_centered_matrix);
    S = MultiplyMatrix(transpose_mean_centered_matrix, mean_centered_matrix);
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            S[i][j] = (E1_EN)*S[i][j];          
        }                                  
    }

    //step3&4  
     for(int h = 0; h < dimension - 1; h++){    
        result = DomEigenvector(S, number_of_iteration - 1);

    	eigenvector2 = MultiplyMatrix(S, result);       
         
        for (int j = 0; j < matrix[0].size(); j++){
            result_temp[j][0] = abs(eigenvector2[j][0]);
        }
        
        for (int i = 0; i < matrix[0].size()-1; i++){
            if(result_temp[i][0] > result_temp[i+1][0]){
                result_temp[i+1][0] = result_temp[i][0];
            }
        }     
	    float max_value = result_temp[matrix[0].size()-1][0];

        float norm_constant = (float)1/max_value;
    
        for (int i = 0; i < matrix[0].size(); i++){
                eigenvector[i][0] = eigenvector2[i][0]*norm_constant;        /////////this is the dominant eigenvector of S. doing one more iteration on result.
        }    

        float num = eigenvector2[0][0] / result[0][0];       ///eigenvalue end///// eigenvalue = num

        float sqsum = 0;
        for(int i = 0; i < matrix[0].size(); i++){
            sqsum += eigenvector[i][0]*eigenvector[i][0];
        }
        float sqroot = sqrt(sqsum);                         // sqareroot of sqsum

        for(int i = 0; i < matrix[0].size(); i++){
            normalized[i][0] = (1/sqroot)*eigenvector[i][0];          
        }
        transpose_normalized = TransposeMatrix(normalized);

        tmp = MultiplyMatrix(normalized, transpose_normalized);
        for(int i = 0; i < matrix[0].size(); i++){
            for(int j = 0; j < matrix[0].size(); j++){
                tmp[i][j] = num*tmp[i][j];
            }
        }

        S = SubtractMatrix(S, tmp);                            ///////deflated matrix

        for(int i = 0; i < matrix[0].size(); i++){
            reduced_basisMatrix[i][h] = eigenvector[i][0];
        }
    }

    result = DomEigenvector(S, number_of_iteration);
    for(int i = 0; i < matrix[0].size(); i++){
        reduced_basisMatrix[i][dimension - 1] = result[i][0];
    }


    cout << "reduced basis matrix is\n";
    for(int i = 0; i < matrix[0].size(); i++){
        for(int j = 0; j < dimension; j++){                           
            cout << reduced_basisMatrix[i][j] << "  ";
        }
        cout << "\n";
    }

    //step5
    plainres = MultiplyMatrix(matrix, reduced_basisMatrix);

    cout << "\n reduced matrix is \n";

    return plainres;
}










void newPCAMatrix_step1(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
	vector<LweSample*> centering(matrix[0].size());
	for(int i = 0; i < matrix[0].size(); i++){
		centering[i] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	}
    LweSample* P2C_N = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* P2C_1 = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* P2C_1_N = new_gate_bootstrapping_ciphertext_array(length, bk->params);

    // Compute the mean of each variable   // μ = (1/n) * ∑(i=1 to n) Xi
    for(int j = 0; j < matrix[0].size(); j++){
        for(int k = 0; k < length; k++){
            bootsCONSTANT(&centering[j][k], 0, bk);
        }
    }
    newP2C(P2C_1, 1, length, bk);
    newP2C(P2C_N, matrix.size(), length, bk);
    newRealDiv(P2C_1_N, P2C_1, P2C_N, length, bk);


    for(int j = 0; j < matrix[0].size(); j++){
        for(int k = 0; k < length; k++){
            bootsCOPY(&centering[j][k], &matrix[0][j][k], bk);
        }
        for(int i = 1; i < matrix.size(); i++){
            newADD(centering[j], matrix[i][j], centering[j], length, bk);
        }
    }
    for(int j = 0; j < matrix[0].size(); j++){
        for(int k = 0; k < length; k++){
            bootsCOPY(&res[0][j][k], &centering[j][k], bk);
        }
    }
    for(int i = 0; i < matrix[0].size(); i++){
        newMultiReal(res[0][i], centering[i], P2C_1_N, length, bk);
    }

    for(int i = 0; i < matrix[0].size(); i++){
		delete_gate_bootstrapping_ciphertext_array(length, centering[i]);
	}
    delete_gate_bootstrapping_ciphertext_array(length, P2C_N);
    delete_gate_bootstrapping_ciphertext_array(length, P2C_1);
    delete_gate_bootstrapping_ciphertext_array(length, P2C_1_N);
}


void newPCAMatrix_step2(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, vector<LweSample*> centering, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    LweSample* P2C_N = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* P2C_1 = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* P2C_1_N = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	vector<vector<LweSample*>> mean_centered_matrix(matrix.size(), vector<LweSample*>(matrix[0].size())); //NxM
	vector<vector<LweSample*>> transpose_mean_centered_matrix(matrix[0].size(), vector<LweSample*>(matrix.size())); //MxN
	for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
		    mean_centered_matrix[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
            transpose_mean_centered_matrix[j][i] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
	}
	vector<vector<LweSample*>> S(matrix[0].size(), vector<LweSample*>(matrix[0].size())); //MxM
	for(int i = 0; i < matrix[0].size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
		    S[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
	}

    // Center the data by subtracting the mean from each observation   // Z = X - μ
    newP2C(P2C_1, 1, length, bk);
    newP2C(P2C_N, matrix.size(), length, bk);
    newRealDiv(P2C_1_N, P2C_1, P2C_N, length, bk);

    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            newSUB(mean_centered_matrix[i][j], matrix[i][j], centering[j], length, bk);
        }
    }

    // Compute the covariance matrix   // S = (1/n) * ZᵀZ
    newTransposeMatrix(transpose_mean_centered_matrix, mean_centered_matrix, length, bk);
    newmultRealMatrix(S, transpose_mean_centered_matrix, mean_centered_matrix, length, bk);
    newScalarMultMatrix2(res, S, P2C_1_N, length, bk);
	
    delete_gate_bootstrapping_ciphertext_array(length, P2C_N);
    delete_gate_bootstrapping_ciphertext_array(length, P2C_1);
    delete_gate_bootstrapping_ciphertext_array(length, P2C_1_N);
	for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
		    delete_gate_bootstrapping_ciphertext_array(length, mean_centered_matrix[i][j]);
		    delete_gate_bootstrapping_ciphertext_array(length, transpose_mean_centered_matrix[j][i]);
        }
	}
	for(int i = 0; i < matrix[0].size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
		    delete_gate_bootstrapping_ciphertext_array(length, S[i][j]);
        }
	}
}


void newPCAMatrix_step4(vector<vector<LweSample*>> res, vector<vector<LweSample*>> eigenvalue_S, vector<vector<LweSample*>> eigenvector_S, const int reduced_basis, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    LweSample* temp = new_gate_bootstrapping_ciphertext_array(2, bk->params);

    // Sort the eigenvalues in descending order and arrange the corresponding eigenvectors accordingly
    // λ₁ ≥ λ₂ ≥ ... ≥ λp, V = [v₁, v₂, ..., vp]
    for (int h = 0; h < eigenvalue_S.size() - 1; h++) {
        for (int j = 0; j < eigenvalue_S.size() - 1 - h; j++) {
            newCompLE(&temp[0], eigenvalue_S[j][0], eigenvalue_S[j+1][0], length, bk); // eigenvalue_S[j] >= eigenvalue_S[j+1] -> temp[0] = 1
            bootsNOT(&temp[1], &temp[0], bk);
            for(int i = 0; i < eigenvector_S.size(); i++){
                for(int k = 0; k < length; k++){
                    bootsMUX(&eigenvector_S[i][j][k], &temp[0], &eigenvector_S[i][j][k], &eigenvector_S[i][j+1][k], bk); // max
                    bootsMUX(&eigenvector_S[i][j+1][k], &temp[1], &eigenvector_S[i][j][k], &eigenvector_S[i][j+1][k], bk); // min
                    bootsMUX(&eigenvalue_S[j][0][k], &temp[0], &eigenvalue_S[j][0][k], &eigenvalue_S[j+1][0][k], bk); // max
                    bootsMUX(&eigenvalue_S[j+1][0][k], &temp[1], &eigenvalue_S[j][0][k], &eigenvalue_S[j+1][0][k], bk); // min
                }
            }
        }
    }

    // W[][] = [v₁, v₂, ..., vr]   // reduced basis
    for(int i = 0; i < eigenvector_S.size(); i++){
        for(int j = 0; j < reduced_basis; j++){
            for(int k = 0; k < length; k++){
                bootsCOPY(&res[i][j][k], &eigenvector_S[i][j][k], bk);
            }
        }
    }

    delete_gate_bootstrapping_ciphertext_array(2, temp);
}


void newPCAMatrix_step5(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, vector<vector<LweSample*>> reduced_Matrix, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    vector<vector<LweSample*>> transpose_reducedMatrix(reduced_Matrix[0].size(), vector<LweSample*>(reduced_Matrix.size()));
	for(int i = 0; i < reduced_Matrix[0].size(); i++){
        for(int j = 0; j < reduced_Matrix.size(); j++){
            transpose_reducedMatrix[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
	}
    vector<vector<LweSample*>> transpose_aMatrix(aMatrix[0].size(), vector<LweSample*>(aMatrix.size()));
	for(int i = 0; i < aMatrix[0].size(); i++){
        for(int j = 0; j < aMatrix.size(); j++){
            transpose_aMatrix[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
	}

    // Project the data onto the new subspace by computing Z' = ZW   // compute reduced dimensionality data
    newTransposeMatrix(transpose_reducedMatrix, reduced_Matrix, length, bk);
    newTransposeMatrix(transpose_aMatrix, aMatrix, length, bk);
    newmultRealMatrix(res, transpose_reducedMatrix, transpose_aMatrix, length, bk);

	for(int i = 0; i < reduced_Matrix[0].size(); i++){
        for(int j = 0; j < reduced_Matrix.size(); j++){
		    delete_gate_bootstrapping_ciphertext_array(length, transpose_reducedMatrix[i][j]);
        }
	}
	for(int i = 0; i < aMatrix[0].size(); i++){
        for(int j = 0; j < aMatrix.size(); j++){
		    delete_gate_bootstrapping_ciphertext_array(length, transpose_aMatrix[i][j]);
        }
	}	
}




void newPCAMatrix_step3_4_reducedbasis(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, const int number_of_iteration, const int dimension, const int length, const TFheGateBootstrappingCloudKeySet* bk) {
    vector<vector<LweSample*>> result(matrix.size(), vector<LweSample*>(1));
    vector<vector<LweSample*>> eigenvector(matrix.size(), vector<LweSample*> (1));
    vector<vector<LweSample*>> eigenvector2(matrix.size(), vector<LweSample*> (1));
    vector<vector<LweSample*>> result_temp(matrix.size(), vector<LweSample*> (1));
    for(int i = 0; i < matrix.size(); i++){
        result[i][0] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        eigenvector[i][0] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        eigenvector2[i][0] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        result_temp[i][0] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    }
    LweSample* max_value = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* norm_constant = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* denominator = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* num = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* P2C_1 = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    vector<vector<LweSample*>> tmp(matrix.size(), vector<LweSample*>(matrix.size()));
	for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix.size(); j++){
            tmp[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
	}
    vector<vector<LweSample*>> reduced_basisMatrix(matrix.size(), vector<LweSample*>(dimension));        
	for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < dimension; j++){
            reduced_basisMatrix[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
	}
    LweSample* sqsum = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* square = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* sqroot = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    vector<vector<LweSample*>> normalized(matrix.size(), vector<LweSample*>(1));
    vector<vector<LweSample*>> transpose_normalized(1, vector<LweSample*>(matrix.size()));
    for(int i = 0; i < matrix.size(); i++){
        normalized[i][0] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        transpose_normalized[0][i] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    }

    newP2C(P2C_1, 1, length, bk);
    for(int h = 0; h < dimension - 1; h++){                                                      
        newDomEigenvector(result, matrix, number_of_iteration - 1, length, bk);
        newmultRealMatrix(eigenvector2, matrix, result, length, bk); 
        for (int j = 0; j < matrix.size(); j++){
            newABS(result_temp[j][0], eigenvector2[j][0], length, bk);
        }
        newMaxValue(max_value, result_temp, length, bk);  
        newRealDiv(norm_constant, P2C_1, max_value, length, bk);
        newScalarMultMatrix2(eigenvector, eigenvector2, norm_constant, length, bk);  /////////this is the dominant eigenvector of S. doing one more iteration on result.

        newRealDiv(denominator, P2C_1, result[0][0], length, bk);
        newMultiReal(num, eigenvector2[0][0], denominator, length, bk);          // eigenvalue = num // eigenvalue associated with 1st dominant eigenvector

        ////deflated matrix -> matrix
        for(int k = 0; k < length; k++){
            bootsCONSTANT(&sqsum[k], 0, bk);
        }
        for(int i = 0; i < matrix.size(); i++){
            newMultiReal(square, eigenvector[i][0], eigenvector[i][0], length, bk);
            newADD(sqsum, sqsum, square, length, bk);
        }
        HomSqroot(sqroot, sqsum, length, bk);                         // sqareroot of sqsum ////////////////////////////////////////////
        newRealDiv(denominator, P2C_1, sqroot, length, bk);
        for(int i = 0; i < matrix.size(); i++){
            newMultiReal(normalized[i][0], denominator, eigenvector[i][0], length, bk);        
        }
        newTransposeMatrix(transpose_normalized, normalized, length, bk);
        newmultRealMatrix(tmp, normalized, transpose_normalized, length, bk);

        #pragma omp parallel for
        for(int i = 0; i < matrix.size(); i++){
            for(int j = 0; j < matrix.size(); j++){
                newMultiReal(tmp[i][j], num, tmp[i][j], length, bk);        
            }
        }

        newSubMatrix(matrix, matrix, tmp, length, bk);                          ///////deflated matrix

        for(int i = 0; i < matrix.size(); i++){
            for(int k = 0; k < length; k++){
                bootsCOPY(&reduced_basisMatrix[i][h][k], &eigenvector[i][0][k], bk);
            }
        }
    }   

    newDomEigenvector(result, matrix, number_of_iteration, length, bk);         /////////this is the eigenvector of S
    for(int i = 0; i < matrix.size(); i++){
        for(int k = 0; k < length; k++){
            bootsCOPY(&reduced_basisMatrix[i][dimension - 1][k], &result[i][0][k], bk);                 
        }
    }

    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < dimension; j++){                                                 
            for(int k = 0; k < length; k++){
                bootsCOPY(&res[i][j][k], &reduced_basisMatrix[i][j][k], bk);                  
            }
        }
    }
    

    for(int i = 0; i < matrix.size(); i++){
        delete_gate_bootstrapping_ciphertext_array(length, result[i][0]);
        delete_gate_bootstrapping_ciphertext_array(length, eigenvector[i][0]);
        delete_gate_bootstrapping_ciphertext_array(length, eigenvector2[i][0]);
        delete_gate_bootstrapping_ciphertext_array(length, result_temp[i][0]);
        delete_gate_bootstrapping_ciphertext_array(length, normalized[i][0]);
        delete_gate_bootstrapping_ciphertext_array(length, transpose_normalized[0][i]);
    }
    delete_gate_bootstrapping_ciphertext_array(length, max_value);
    delete_gate_bootstrapping_ciphertext_array(length, norm_constant);
    delete_gate_bootstrapping_ciphertext_array(length, denominator);
    delete_gate_bootstrapping_ciphertext_array(length, num);
    delete_gate_bootstrapping_ciphertext_array(length, P2C_1);
    delete_gate_bootstrapping_ciphertext_array(length, sqsum);
    delete_gate_bootstrapping_ciphertext_array(length, square);
    delete_gate_bootstrapping_ciphertext_array(length, sqroot);
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix.size(); j++){
            delete_gate_bootstrapping_ciphertext_array(length, tmp[i][j]);
        }
        for(int j = 0; j < dimension; j++){                                                      
		    delete_gate_bootstrapping_ciphertext_array(length, reduced_basisMatrix[i][j]);
        }
    }
}





void newPCAMatrix(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, const int number_of_iteration, const int dimension, const int length, const TFheGateBootstrappingCloudKeySet* bk){
	vector<LweSample*> centering(matrix[0].size());
	for(int i = 0; i < matrix[0].size(); i++){
		centering[i] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
	}
    LweSample* P2C_N = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* P2C_1 = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* P2C_1_N = new_gate_bootstrapping_ciphertext_array(length, bk->params);

	vector<vector<LweSample*>> mean_centered_matrix(matrix.size(), vector<LweSample*>(matrix[0].size())); //NxM
	vector<vector<LweSample*>> transpose_mean_centered_matrix(matrix[0].size(), vector<LweSample*>(matrix.size())); //MxN
	for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
		    mean_centered_matrix[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
            transpose_mean_centered_matrix[j][i] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
	}
	vector<vector<LweSample*>> S(matrix[0].size(), vector<LweSample*>(matrix[0].size())); //MxM
	for(int i = 0; i < matrix[0].size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
		    S[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
	}
    vector<vector<LweSample*>> tmp(matrix[0].size(), vector<LweSample*>(matrix[0].size()));
	for(int i = 0; i < matrix[0].size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            tmp[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
	}
    vector<vector<LweSample*>> result(matrix[0].size(), vector<LweSample*>(1));
    vector<vector<LweSample*>> eigenvector(matrix[0].size(), vector<LweSample*>(1));
    vector<vector<LweSample*>> eigenvector2(matrix[0].size(), vector<LweSample*>(1));
    vector<vector<LweSample*>> result_temp(matrix[0].size(), vector<LweSample*>(1));
	for(int i = 0; i < matrix[0].size(); i++){
        for(int j = 0; j < 1; j++){
            result[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
            eigenvector[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
            eigenvector2[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
            result_temp[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
	}
    vector<vector<LweSample*>> reduced_BasisMatrix(matrix[0].size(), vector<LweSample*>(dimension));          
	for(int i = 0; i < matrix[0].size(); i++){
        for(int j = 0; j < dimension; j++){
            reduced_BasisMatrix[i][j] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        }
	}
    LweSample* max_value = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* norm_constant = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* denominator = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* num = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* sqsum = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* square = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    LweSample* sqroot = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    vector<vector<LweSample*>> normalized(matrix[0].size(), vector<LweSample*>(1));
    vector<vector<LweSample*>> transpose_normalized(1, vector<LweSample*>(matrix[0].size()));
    for(int i = 0; i < matrix[0].size(); i++){
        normalized[i][0] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
        transpose_normalized[0][i] = new_gate_bootstrapping_ciphertext_array(length, bk->params);
    }

    // step1
    // Compute the mean of each variable   // μ = (1/n) * ∑(i=1 to n) Xi
    std::chrono::system_clock::time_point start1 = std::chrono::system_clock::now();
    for(int j = 0; j < matrix[0].size(); j++){
        for(int k = 0; k < length; k++){
            bootsCOPY(&centering[j][k], &matrix[0][j][k], bk);
        }
        for(int i = 1; i < matrix.size(); i++){
            newADD(centering[j], matrix[i][j], centering[j], length, bk);
        }
    }
    newP2C(P2C_1, 1, length, bk);
    newP2C(P2C_N, matrix.size(), length, bk);
    newRealDiv(P2C_1_N, P2C_1, P2C_N, length, bk);
    for(int i = 0; i < matrix[0].size(); i++){
        newMultiReal(centering[i], centering[i], P2C_1_N, length, bk);
    }
    std::chrono::duration<float> sec1 = std::chrono::system_clock::now() - start1;
    std::cout << "(in Chrono) step1 done in " << sec1.count() << " seconds...\n" << std::endl;

    // step2
    // Center the data by subtracting the mean from each observation   // Z = X - μ
    std::chrono::system_clock::time_point start2 = std::chrono::system_clock::now();
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            newSUB(mean_centered_matrix[i][j], matrix[i][j], centering[j], length, bk);
        }
    }
    // Compute the covariance matrix   // S = (1/n) * ZᵀZ
    newTransposeMatrix(transpose_mean_centered_matrix, mean_centered_matrix, length, bk);
    newmultRealMatrix(S, transpose_mean_centered_matrix, mean_centered_matrix, length, bk);
    newScalarMultMatrix2(S, S, P2C_1_N, length, bk);
    std::chrono::duration<float> sec2 = std::chrono::system_clock::now() - start2;
    std::cout << "(in Chrono) step2 done in " << sec2.count() << " seconds...\n" << std::endl;

    //step3&4    
    std::chrono::system_clock::time_point start3 = std::chrono::system_clock::now();
    for(int h = 0; h < dimension - 1; h++){                     
        std::chrono::system_clock::time_point start3_1 = std::chrono::system_clock::now();                                      
        newDomEigenvector(result, S, number_of_iteration - 1, length, bk);
        newmultRealMatrix(eigenvector2, S, result, length, bk);  
        for (int j = 0; j < matrix[0].size(); j++){
            newABS(result_temp[j][0], eigenvector2[j][0], length, bk);
        }
        newMaxValue(max_value, result_temp, length, bk);  
        newRealDiv(norm_constant, P2C_1, max_value, length, bk);
        newScalarMultMatrix2(eigenvector, eigenvector2, norm_constant, length, bk);  /////////this is the dominant eigenvector of S. doing one more iteration on result.
        std::chrono::duration<float> sec3_1 = std::chrono::system_clock::now() - start3_1;
        std::cout << "(in Chrono) step3_1 done in " << sec3_1.count() << " seconds...\n" << std::endl;
        
        std::chrono::system_clock::time_point start3_2 = std::chrono::system_clock::now();  
        newRealDiv(num, eigenvector2[0][0], result[0][0], length, bk);    // eigenvalue = num // eigenvalue associated with 1st dominant eigenvector

        for(int k = 0; k < length; k++){
            bootsCONSTANT(&sqsum[k], 0, bk);
        }
        for(int i = 0; i < matrix[0].size(); i++){
            newMultiReal(square, eigenvector[i][0], eigenvector[i][0], length, bk);
            newADD(sqsum, sqsum, square, length, bk);
        }

        HomSqroot(sqroot, sqsum, length, bk);                         // sqareroot of sqsum ////////////////////////////////////////////

        newRealDiv(denominator, P2C_1, sqroot, length, bk);

        for(int i = 0; i < matrix[0].size(); i++){
            newMultiReal(normalized[i][0], denominator, eigenvector[i][0], length, bk);        
        }
        std::chrono::duration<float> sec3_2 = std::chrono::system_clock::now() - start3_2;
        std::cout << "(in Chrono) step3_2 done in " << sec3_2.count() << " seconds...\n" << std::endl;

        std::chrono::system_clock::time_point start3_3 = std::chrono::system_clock::now(); 
        newTransposeMatrix(transpose_normalized, normalized, length, bk);

        newmultRealMatrix(tmp, normalized, transpose_normalized, length, bk);

        #pragma omp parallel for
        for(int i = 0; i < matrix[0].size(); i++){
            for(int j = 0; j < matrix[0].size(); j++){
                newMultiReal(tmp[i][j], num, tmp[i][j], length, bk);        
            }
        }

        newSubMatrix(S, S, tmp, length, bk);                          ///////deflated matrix

        for(int i = 0; i < matrix[0].size(); i++){
            for(int k = 0; k < length; k++){
                bootsCOPY(&reduced_BasisMatrix[i][h][k], &eigenvector[i][0][k], bk);
            }
        }
        std::chrono::duration<float> sec3_3 = std::chrono::system_clock::now() - start3_3;
        std::cout << "(in Chrono) step3_3 done in " << sec3_3.count() << " seconds...\n" << std::endl;
    }
    std::chrono::duration<float> sec3 = std::chrono::system_clock::now() - start3;
    std::cout << "(in Chrono) step3 done in " << sec3.count() << " seconds...\n" << std::endl;

    std::chrono::system_clock::time_point start4 = std::chrono::system_clock::now();
    newDomEigenvector(result, S, number_of_iteration, length, bk);         /////////this is the eigenvector of S
    for(int i = 0; i < matrix[0].size(); i++){
        for(int k = 0; k < length; k++){
            bootsCOPY(&reduced_BasisMatrix[i][dimension - 1][k], &result[i][0][k], bk);                 
        }
    }
    std::chrono::duration<float> sec4 = std::chrono::system_clock::now() - start4;
    std::cout << "(in Chrono) step4 done in " << sec4.count() << " seconds...\n" << std::endl;

    //step5
    // Project the data onto the new subspace by computing Z' = ZW   // compute reduced dimensionality data
    std::chrono::system_clock::time_point start5 = std::chrono::system_clock::now();
    newmultRealMatrix(res, matrix, reduced_BasisMatrix, length, bk);
    std::chrono::duration<float> sec5 = std::chrono::system_clock::now() - start5;
    std::cout << "(in Chrono) step5 done in " << sec5.count() << " seconds...\n" << std::endl;


    delete_gate_bootstrapping_ciphertext_array(length, P2C_N);
    delete_gate_bootstrapping_ciphertext_array(length, P2C_1);
    delete_gate_bootstrapping_ciphertext_array(length, P2C_1_N);
	for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
		    delete_gate_bootstrapping_ciphertext_array(length, mean_centered_matrix[i][j]);
		    delete_gate_bootstrapping_ciphertext_array(length, transpose_mean_centered_matrix[j][i]);
        }
	}
    for(int i = 0; i < matrix[0].size(); i++){
		delete_gate_bootstrapping_ciphertext_array(length, centering[i]);
        for(int j = 0; j < matrix[0].size(); j++){
		    delete_gate_bootstrapping_ciphertext_array(length, S[i][j]);
        }
        for(int j = 0; j < 1; j++){
		    delete_gate_bootstrapping_ciphertext_array(length, result[i][j]);
            delete_gate_bootstrapping_ciphertext_array(length, eigenvector[i][j]);
            delete_gate_bootstrapping_ciphertext_array(length, eigenvector2[i][j]);
            delete_gate_bootstrapping_ciphertext_array(length, result_temp[i][j]);
        }
        for(int j = 0; j < dimension; j++){                                                         
		    delete_gate_bootstrapping_ciphertext_array(length, reduced_BasisMatrix[i][j]);
        }
	}
    delete_gate_bootstrapping_ciphertext_array(length, max_value);
    delete_gate_bootstrapping_ciphertext_array(length, norm_constant);
    delete_gate_bootstrapping_ciphertext_array(length, denominator);
    delete_gate_bootstrapping_ciphertext_array(length, num);
	for(int i = 0; i < matrix[0].size(); i++){                                                  
        for(int j = 0; j < matrix[0].size(); j++){
		    delete_gate_bootstrapping_ciphertext_array(length, tmp[i][j]);
        }
	}
    delete_gate_bootstrapping_ciphertext_array(length, sqsum);
    delete_gate_bootstrapping_ciphertext_array(length, square);
    delete_gate_bootstrapping_ciphertext_array(length, sqroot);
    for(int i = 0; i < matrix[0].size(); i++){
        delete_gate_bootstrapping_ciphertext_array(length, normalized[i][0]);
        delete_gate_bootstrapping_ciphertext_array(length, transpose_normalized[0][i]);
    }
}