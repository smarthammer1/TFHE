g++ testmatrixstep1.cpp -o testmatrixstep1 -lm -ltfhe-spqlios-fma -fopenmp
g++ testmatrixstep2_matrix.cpp -o testmatrixstep2_matrix -lm -ltfhe-spqlios-fma -fopenmp
g++ testmatrixstep3_matrix.cpp -o testmatrixstep3_matrix -lm -ltfhe-spqlios-fma -fopenmp
  
chmod 755 testmatrixstep1
chmod 755 testmatrixstep2_matrix
chmod 755 testmatrixstep3_matrix
  
./testmatrixstep1
./testmatrixstep2_matrix
./testmatrixstep3_matrix
