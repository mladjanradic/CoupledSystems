clear, close, clc;

FN = '/Users/mladi/Desktop/GAMMA/Std_Err_Est_part/';

FN_save = '/Users/mladi/Desktop/GAMMA/Hier_Err_Est_part/';

mkdir(FN_save);

i=1;
data = load([FN 'truth_err' num2str(i) '.mat']);
m=data.m;
model_data = data.d;
M_train = data.M_train;

K1 = blkdiag(data.d.LGG2,0*data.d.LGG2);
K2 = blkdiag(0*data.d.LGG2,data.d.LGG2);

clear data;
D=3;

for i=1:15-D
    dataN = load([FN 'truth_err' num2str(i) '.mat']);
    red_dataN = dataN.red_data;
    
    dataM = load([FN 'truth_err' num2str(i+D) '.mat']);
    red_dataM = dataM.red_data;
    
    zaehler = dataM.truth_err_OMEGA1;
    nenner = dataN.truth_err_OMEGA1;
    
    nenner(dataN.ind_seq) = [];
    zaehler(dataN.ind_seq) = [];
    
    ind = find(nenner>(norm(nenner)/100));
    MAX = max(zaehler(ind)./nenner(ind));
    theta(i,1) = MAX(1);
    
    
    %OMEGA2:
    zaehler = dataM.truth_err_OMEGA2;
    nenner = dataN.truth_err_OMEGA2;
    
    nenner(dataN.ind_seq) = [];
    zaehler(dataN.ind_seq) = [];
    
    ind = find(nenner>(norm(nenner)/100));
    MAX = max(zaehler(ind)./nenner(ind));
    theta(i,2) = MAX(1);
    
end

theta



%%


for i=1:15-D
    dataN = load([FN 'truth_err' num2str(i) '.mat']);
    red_dataN = dataN.red_data;
    
    dataM = load([FN 'truth_err' num2str(i+D) '.mat']);
    red_dataM = dataM.red_data;
    d = dataN.d;
    
    M=500;
    M_test = rand_uniform(500, m.model.mu_ranges);
    M_test(1,:) = 1;
    M_test(2,:) = 3;
    
    for k = 1:M
        mu = M_test(:,k);
        m.model = m.model.set_mu(m.model,mu);
        
        m.operators = @operators;
        red_simN = stabilized_lin_stat_rb_simulation(m,red_dataN);
        uN = red_simN.uN;
        red_simM = stabilized_lin_stat_rb_simulation(m,red_dataM);
        uM = red_simN.uN;
        
        UN = dataN.d.RB*red_simN.uN;
        UM = dataM.d.RB*red_simM.uN;
        
        diff = UN - UM;
        DELTA1(i,k) = sqrt(abs(diff'*K1*diff));
        DELTA2(i,k) = sqrt(abs(diff'*K2*diff));
        
        [tru1,tru2,u1,u2] = compute_trace_solutions(m.model,d,d.utilde1,d.utilde2);
        uh = [tru1;tru2];
        U = UN;
        diff = uh - U;
        
        truth_err_OMEGA1(i,k) = sqrt(abs(diff'*K1*diff));
        truth_err_OMEGA2(i,k) = sqrt(abs(diff'*K2*diff));
        
        
    end
    hier_err1(i) = sum(DELTA1(i,:)/(1-theta(i,1)))/M;
    hier_err2(i) = sum(DELTA2(i,:)/(1-theta(i,2)))/M;
    
    ERR1(i) = sum(truth_err_OMEGA1(i,:))/M;
    ERR2(i) = sum(truth_err_OMEGA2(i,:))/M;
    
    save([FN_save,'data' num2str(i) '.mat']);
    disp(i);
end


header = {'N','Hier1','Hier2','ErrOMEGA1','ErrOMEGA2'};
matrix(1:15-D,1) = 1:15-D;
matrix(1:15-D,2) = hier_err1(:);
matrix(1:15-D,3) = hier_err2(:);
matrix(1:15-D,4) = ERR1(:);
matrix(1:15-D,5) = ERR2(:);
matrix2txt(header, matrix, [FN_save 'hier_vs_error.txt'])
