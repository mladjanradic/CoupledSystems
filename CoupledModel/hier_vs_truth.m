clear, close, clc;

D=3;
M=1000;

for i=1:15-D
    data = load(['/Users/mladi/Desktop/HighDim/HierErr_D' num2str(D) ...
        '/data' num2str(i) '.mat']);
    hier1(i,1) = sum(data.hier_err_est1)/M;
    hier2(i,1) = sum(data.hier_err_est2)/M;
    
    err1(i,1) = sum(data.truth_err1)/M;
    err2(i,1) = sum(data.truth_err2)/M;
end

FN = ['/Users/mladi/Desktop/HighDim/HierErr_D' num2str(D) '/'];
header = {'N','HIER1','HIER2','ErrOMEGA1','ErrOMEGA2'};
matrix(1:15-D,1) = 1:15-D;
matrix(1:15-D,2) = hier1;
matrix(1:15-D,3) = hier2;
matrix(1:15-D,4) = err1;
matrix(1:15-D,5) = err2;
matrix
matrix2txt(header, matrix, [FN 'hier_vs_error.txt'])


