#include "newfile.cpp"
#include "newfilepca.cpp"
#include <omp.h>
#include <chrono>

////////////////////////////////////////
/*vector<vector<float>> res is matrix*/

int main() {
    printf("============ Step2 start ============\n");
    //reads the cloud key from file
    FILE* cloud_key = fopen("cloud.key","rb");
    TFheGateBootstrappingCloudKeySet* bk = new_tfheGateBootstrappingCloudKeySet_fromFile(cloud_key);
    fclose(cloud_key);
    //if necessary, the params are inside the key
    const TFheGateBootstrappingParameterSet* params = bk->params;

    //parameter
    int length = 32;
    int res_length = 32;
    // const int d = 2;
    // const int r = 1;
    const int number_of_iteration = 3;      // at least 3
    float num = 3.6;
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

    //room for the ciphertexts
    int row1 = aMatrix.size();
    int col1 = aMatrix[0].size();
    vector<vector<LweSample*>> cipher1(row1, vector<LweSample*>(col1)); 
    // LweSample* cipher1[10][10];
    for (int i = 0; i < aMatrix.size(); i++){
        for (int j = 0; j < aMatrix[0].size(); j++){
            cipher1[i][j] = new_gate_bootstrapping_ciphertext_array(length, params);
        }
    }

    int row2 = bMatrix.size();
    int col2 = bMatrix[0].size();
    vector<vector<LweSample*>> cipher2(row2, vector<LweSample*>(col2)); 
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

    ////////////////res matrix
    vector<vector<LweSample*>> res(row1, vector<LweSample*>(dimension));            ////////////////////        answer      //////////////////////
    for (int i = 0; i < res.size(); i++){
        for (int j = 0; j < res[0].size(); j++){
            res[i][j] = new_gate_bootstrapping_ciphertext_array(length, params);
        }
    }


    //reads the ciphertexts from the cloud file 
    FILE* cloud_data = fopen("cloud.data","rb");


    for (int i = 0; i < aMatrix.size(); i++){
        for (int j = 0; j < aMatrix[0].size(); j++){
            for (int k = 0; k < length; k++){
                import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &cipher1[i][j][k], params);
            }
        }
    }
    for (int i = 0; i < bMatrix.size(); i++){
        for (int j = 0; j < bMatrix[0].size(); j++){
            for (int k = 0; k < length; k++){
                import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &cipher2[i][j][k], params);
            }
        }
    }
    for (int i = 0; i < cMatrix.size(); i++){
        for (int j = 0; j < cMatrix[0].size(); j++){
            for (int k = 0; k < length; k++){
                import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &cipher3[i][j][k], params);
            }
        }
    }
    for (int i = 0; i < dMatrix.size(); i++){
        for (int j = 0; j < dMatrix[0].size(); j++){
            for (int k = 0; k < length; k++){
                import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &cipher4[i][j][k], params);
            }
        }
    }
    for (int i = 0; i < centering.size(); i++){
        for (int k = 0; k < length; k++){
            import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &cipher5[i][k], params);
        }
    }
    for (int k = 0; k < length; k++){
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &cipher6[k], params);
    }

    fclose(cloud_data);

    float time = -clock();
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();


                              
    //////// HE operation
    
    // cout << "### Output Data ###" << endl;
    // cout << "=== Matrix _ PCA === \n"; 

    ///////////      print plainres     ///////////////
    vector<vector<float>> plainres(row1, vector<float>(dimension));              ////////////////////        answer      //////////////////////
    // vector<float> plainres(col1);      
    //////////////////////////

    // plainres = print_matrix(bMatrix);
    // plainres = AddMatrix(aMatrix, bMatrix);
    // plainres = SubtractMatrix(aMatrix, bMatrix);
    // plainres = MultiplyMatrix(aMatrix, bMatrix);
    // plainres = TransposeMatrix(bMatrix);
    // // InverseMatrix(aMatrix, number_of_iteration);
    // plainres = DomEigenvector(aMatrix, number_of_iteration);


    // plainres = PCAMatrix_step1(aMatrix);
    // plainres = PCAMatrix_step2(aMatrix, centering);
    // plainres = PCAMatrix_step4(cMatrix,bMatrix, dimension);
    // plainres = PCAMatrix_step5(bMatrix, dMatrix);
    // plainres = PCAMatrix_step3_4_reducedbasis(aMatrix, number_of_iteration, dimension);

    plainres = PCAMatrix(aMatrix, number_of_iteration, dimension);

    ///////////      print plainres     ///////////////
    for (int i = 0; i < plainres.size(); i++){                   
        for (int j = 0; j < plainres[0].size(); j++){
            cout << plainres[i][j] << "  ";
        }
        cout << "\n";
    }
    // for (int i = 0; i < plainres.size(); i++){                   
    //     cout << plainres[i] << "  ";
    // }
    //////////////////////////


    // newprint_Matrix(res, cipher1, length, bk);
    // newAddMatrix(res, cipher2, cipher3, length, bk);
    // newSubMatrix(res, cipher1, cipher2, length, bk);
    // newmultRealMatrix(res, cipher1, cipher2, length, bk);
    // newTransposeMatrix(res, cipher1, length, bk);
    // newDomEigenvector(res, cipher1, number_of_iteration, length, bk);


    // newInverseMatrix(res, cipher1, number_of_iteration, length, bk); //if a=NxM, then res = MxN
    // newP2C(res, num, length, bk);
    // newScalarMultMatrix(res, cipher1, num, length, bk);
    // newScalarMultMatrix2(res, cipher1, cipher5, length, bk);


    // newPCAMatrix_step1(res, cipher1, length, bk);
    // newPCAMatrix_step2(res, cipher1, cipher5, length, bk);
    // newPCAMatrix_step4(res, cipher3, cipher2, dimension, length, bk);
    // newPCAMatrix_step5(res, cipher2, cipher4, length, bk);
    // newPCAMatrix_step3_4_reducedbasis(res, cipher1, number_of_iteration, dimension, length, bk);

    // newPCAMatrix(res, cipher1, number_of_iteration, dimension, length, bk);

    // // newPCAMatrix_step3



    //////// 
    std::chrono::duration<float>sec = std::chrono::system_clock::now() - start;
    time += clock();
    time = time/(CLOCKS_PER_SEC);
    printf("done in %f seconds...\n", time);
    std::cout << "(in Chrono) done in " << sec.count() << " seconds...\n" << std::endl;

    //export ciphertext to a file (for the cloud)
    FILE* answer_data = fopen("answer.data","wb");

    for (int i = 0; i < res.size(); i++){
        for (int j = 0; j < res[0].size(); j++){
            for (int k = 0; k < res_length; k++){ 
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &res[i][j][k], params);
            }
        }
    }

    fclose(answer_data);

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

    for (int i = 0; i < res.size(); i++){
        for (int j = 0; j < res[0].size(); j++){
            delete_gate_bootstrapping_ciphertext_array(length, res[i][j]);
        }
    }
    delete_gate_bootstrapping_cloud_keyset(bk);

    return 0;
}


