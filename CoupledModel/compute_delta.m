function delta = compute_delta(m,d,model_data)
AII1 = d.AII1{1};
AIG1 = d.AIG1{1};
AGI1 = d.AGI1{1};
AGG1 = d.AGG1{1};
SIGMA1 = AGI1*(AII1\AIG1);

AII2 = d.AII2{1};
AIG2 = d.AIG2{1};
AGI2 = d.AGI2{1};
AGG2 = d.AGG2{1};

SIGMA2 = AGI2*(AII2\AIG2);
LGG2 = d.LGG2;

FN='/Users/mladi/Desktop/GAMMA/deltas/';
mkdir(FN);

for i=1:size(d.M_train,2)

mu = d.M_train(:,i);

AGG1 = AGG1+mu(3)*d.LGG1;
AGG2 = AGG2+mu(4)*LGG2;

S1 = AGG1 - SIGMA1;
S2 = AGG2 - SIGMA2;


alpha1 = eigs(0.5*(S1+S1'),d.LGG2,1,'sm');
alpha2 = eigs(0.5*(S2+S2'),LGG2,1,'sm');

gamma1 = eigs(0.5*(S1+S1'),d.LGG2,1,'lm');
gamma2 = eigs(0.5*(S2+S2'),LGG2,1,'lm');

beta1 = mu(3);
beta2 = mu(4);
gammab1 = mu(3);
gammab2 = mu(4);
gammab2b1 = abs(mu(3) - mu(4));

zaehler = (1./alpha2)*(1./(gammab2^2))*gammab2b1;
nenner = alpha2*gamma2*(beta2^2)+alpha1;
delta = zaehler/nenner

disp(i)

save([FN 'delta' num2str(i) '.mat']);

end

end
