function [u1,u2,tru1,tru2] = detailed_simulation_test(model,model_data,d)

mu = get_mu(model);
mu1 = [model.base_model.varepsilon;1;1];
mu2 = [mu(1);mu(2)];
AII1 = 0;AIG1 = 0;AGI1 = 0;AGG1 = 0;
for i=1:size(d.AII1,2)
    AII1 = AII1+mu1(i)*d.AII1{i};
    AIG1 = AIG1+mu1(i)*d.AIG1{i};
    AGI1 = AGI1+mu1(i)*d.AGI1{i};
    AGG1 = AGG1+mu1(i)*d.AGG1{i};
end
AII2=0;AIG2=0;AGI2=0;AGG2=0;
for i=1:size(d.AII2,2)
    AII2 = AII2+ mu2(i)*d.AII2{i};
    AIG2 = AIG2+mu2(i)*d.AIG2{i};
    AGI2 = AGI2+mu2(i)*d.AGI2{i};
    AGG2 = AGG2+mu2(i)*d.AGG2{i};
end
r1=0;rG1=0;r2=0;rG2=0;
for i=1:size(d.r1,2)
    r1 = r1+d.r1{i};
    rG1 = rG1+d.rG1{i};
end
for i=1:size(d.r2,2)
    r2 = r2+d.r2{i};
    rG2 = rG2+d.rG2{i};
end
LGG1 = d.LGG1;
LGG2 = d.LGG2;

utilde1 = AII1\r1;
utilde2 = AII2\r2;
AGG1 = AGG1 - AGI1*(AII1\AIG1);
AGG2 = AGG2 - AGI2*(AII2\AIG2);
rG1 = rG1 - AGI1*utilde1;
rG2 = rG2 - AGI2*utilde2;

AGG = [AGG1+mu(3)*LGG1 -mu(3)*LGG1 ; mu(4)*LGG2 AGG2-mu(4)*LGG2];
r = [rG1;rG2];
%uG can be seen as "[uG1;uG2]"
uG = AGG \ r;

tru1 = uG(1:16);
uG(1:16) = [];
tru2 = uG;

u1 = utilde1 - AII1\(AIG1*tru1);
u2 = utilde2 - AII2\(AIG2*tru2);

