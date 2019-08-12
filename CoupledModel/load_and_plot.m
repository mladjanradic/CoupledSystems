clc
SEQUENCE = zeros(15,1);
for i =1:15
data = load(['/Users/mladi/Desktop/HighDim/Error_Std_Err_Est_3/truth_err' num2str(i) '.mat']);
SEQUENCE(i) = max(data.truth_err);
end
semilogy(SEQUENCE,'r'); hold on; semilogy(data.err_seq,'b');