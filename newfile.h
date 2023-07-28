void printMatrix(vector<vector<float>> matrix);
vector<vector<float>> print_matrix(vector<vector<float>> matrix);
vector<vector<float>> AddMatrix(vector<vector<float>> aMatrix, vector<vector<float>> bMatrix);
vector<vector<float>> SubtractMatrix(vector<vector<float>> aMatrix, vector<vector<float>> bMatrix);
vector<vector<float>> MultiplyMatrix(vector<vector<float>> aMatrix, vector<vector<float>> bMatrix);
vector<vector<float>> TransposeMatrix(vector<vector<float>> matrix);
// vector<vector<float>> InverseMatrix(vector<vector<float>> matrix, const int number_of_iteration);

vector<vector<float>> DomEigenvector(vector<vector<float>> matrix, const int number_of_iteration);
vector<vector<float>> DomEigenvector_QR(vector<vector<float>> matrix, const int number_of_iteration);



void decrypting_code(vector<vector<LweSample*>> res, const int length);
void practice(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, vector<vector<LweSample*>> bMatrix, const int length, const TFheGateBootstrappingCloudKeySet* bk);

void newprint_Matrix(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newAddMatrix(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, vector<vector<LweSample*>> bMatrix, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newmultRealMatrix(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, vector<vector<LweSample*>> bMatrix, const int length, const TFheGateBootstrappingCloudKeySet* bk);
// void newmultRealMatrix2(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, vector<vector<LweSample*>> bMatrix, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newSubMatrix(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, vector<vector<LweSample*>> bMatrix, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newTraceMatrix(LweSample* res, vector<vector<LweSample*>> aMatrix, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newTransposeMatrix(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newP2C(LweSample* res, const float num, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newScalarMultMatrix(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, const float num, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newScalarMultMatrix2(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, const LweSample* b, const int length, const TFheGateBootstrappingCloudKeySet* bk);
// void newInverseMatrix(vector<vector<LweSample*>> res, vector<vector<LweSample*>> aMatrix, const int number_of_iteration, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newDomEigenvector(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, const int number_of_iteration, const int length, const TFheGateBootstrappingCloudKeySet* bk);

//  newEigenvectorMatrix
//  newEigenvalueMatrix
// void newDomEigenvalue(vector<LweSample*> res, vector<vector<LweSample*>> matrix, const int number_of_iteration, const int length, const TFheGateBootstrappingCloudKeySet* bk);
// void newDomEigenvalue(LweSample* res, vector<vector<LweSample*>> matrix, const int number_of_iteration, const int length, const TFheGateBootstrappingCloudKeySet* bk);
void newDomEigenvector_QR(vector<vector<LweSample*>> res, vector<vector<LweSample*>> matrix, const int number_of_iteration, const int length, const TFheGateBootstrappingCloudKeySet* bk);

