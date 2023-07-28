vector<vector<float>> PCAMatrix_step1(vector<vector<float>> matrix);
vector<vector<float>> PCAMatrix_step2(vector<vector<float>> matrix, vector<float> centering);
vector<vector<float>> PCAMatrix_step4(vector<vector<float>> eigenvalue_S, vector<vector<float>> eigenvector_S, const int reduced_basis);
vector<vector<float>> PCAMatrix_step5(vector<vector<float>> aMatrix, vector<vector<float>> reduced_Matrix);
vector<vector<float>> PCAMatrix_step3_4_reducedbasis(vector<vector<float>> matrix, const int number_of_iteration, const int dimension);

vector<vector<float>> PCAMatrix(vector<vector<float>> matrix, const int number_of_iteration, const int dimension);

void newPCAMatrix_step1(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newPCAMatrix_step2(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, vector<LweSample*> centering, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newPCAMatrix_step4(vector<vector<LweSample*>> res, vector<vector<LweSample*>> eigenvalue_S, vector<vector<LweSample*>> eigenvector_S, const int reduced_basis, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newPCAMatrix_step5(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, vector<vector<LweSample*>> reduced_Matrix, const int length, const TFheGateBootstrappingCloudKeySet* bk);

// void newPCAMatrix_step3_4_reducedbasis(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, const int number_of_iteration, const int dimension, const int length, const TFheGateBootstrappingCloudKeySet* bk);
 
void decrypting_code(vector<vector<LweSample*>> res, const int length);

// newPCAMatrix_step3




