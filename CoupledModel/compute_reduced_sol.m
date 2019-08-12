function [uN1,uNG1,uN2,uNG2] = compute_reduced_sol(m,r)

mu = m.mus;
mu1 = [m.varepsilon;1;1];
mu2 = [mu(1);mu(2)];

ANII1=0;ANGI1=0;ANIG1=0;ANGG1=0;
for i=1:size(r.ANII1,2)
    ANII1 = ANII1+mu1(i)*r.ANII1{i};
    ANGI1 = ANGI1+mu1(i)*r.ANGI1{i};
    ANIG1 = ANIG1+mu1(i)*r.ANIG1{i};
    ANGG1 = ANGG1+mu1(i)*r.ANGG1{i};
end
ANII2=0;ANGI2=0;ANIG2=0;ANGG2=0;
for i=1:size(r.ANII2,2)
    ANII2 = ANII2+mu2(i)*r.ANII2{i};
    ANGI2 = ANGI2+mu2(i)*r.ANGI2{i};
    ANIG2 = ANIG2+mu2(i)*r.ANIG2{i};
    ANGG2 = ANGG2+mu2(i)*r.ANGG2{i};
end

rN1 = 0; rNG1=0; rN2=0;rNG2=0;
for i=1:size(r.rN1,2)
    rN1 = rN1 + r.rN1{i};
    rNG1 = rNG1 + r.rNG1{i};
end
for i=1:size(r.rN2,2)
    rN2 = rN2 + r.rN2{i};
    rNG2 = rNG2 + r.rNG2{i};
end

hatu1 = ANII1\rN1;
FG1 = rNG1 - ANGI1*hatu1;
AN1 = ANGG1 - ANGI1*(ANII1\ANIG1);

hatu2 = ANII2\rN2;
FG2 = rNG2 - ANGI2*hatu2;
AN2 = ANGG2 - ANGI2*(ANII2\ANIG2);

A = [AN1+mu(3)*r.LNGG1 -mu(3)*r.LNGG1_cross ; mu(4)*r.LNGG2_cross AN2-mu(4)*r.LNGG2];


N1 = length(rNG1);
N2 = length(rNG2);
uG = A\[FG1;FG2];
uNG1 = uG(1:N1);
uG(1:N1)=[];
uNG2 = uG;

uN1 = hatu1 - ANII1\(ANIG1*uNG1);
uN2 = hatu2 - ANII2\(ANIG2*uNG2);



