function [MT_COEF,HT_COEF,rhoc,cgmixc] = MT_HT_COEF(Tchh)

%% Everything we need for this function:
Molwt=[31.99886, 18.01528, 2.015894, 44.00964, 28.01021, 28.01344];
nu=[ 16.6 , 12.7 , 7.07 , 26.9,  18.9,  17.9 ];
Nsc=5;
sub_thick = 150.0e-06;
thick_L2 = 7.2e-5;
Dch =  1.e-3;
dh = Dch - ((2*sub_thick) + (2*thick_L2) );
hyd_dia = dh;
P = 1.0d5;

%% Declaration & Allocation
Df = zeros(Nsc,1);
MT_COEF = zeros(Nsc,1);

%% Computation:
for iu = 1:Nsc
    nnu = (  (1/Molwt(iu))  + (1/Molwt(Nsc+1)) )^0.5;
    ddn = ( (nu(iu)^(0.333)) + (nu(Nsc+1)^(0.333)) )^2;
    Df(iu) = 0.1*(1.013d-02)*(Tchh^1.75)*(nnu/ddn)*(1/P);
    MT_COEF(iu) = Df(iu)*2.47/hyd_dia;
end

rhoc= 0.8;
cgmixc=1000;
lambda_g = (2.66d-04)*(Tchh^0.805);
HT_COEF = 2.47*lambda_g/hyd_dia;

end





% subroutine MT_HT_COEF(cgmix,rho,MT_COEF,HT_COEf, i_m_f1,Tchh)
%  use SHARED_DATA
% implicit none
% integer ip,jp,kp,iu
% double precision X(Nsc+1),i_m_f(Nsc+1),i_m_f1(Nsc),N2molf
% double precision Cp(Nsc+1)
% double precision dent,Nus_num(Nc),rrh,NM(Nsc+1),sum
% double precision lambda_g, axial_dist(Nc),Tchh
% double precision kme(Nsc*Nc),kmi(Nsc*Nc),Sher_int(Nsc)
% double precision nnu,ddn,Df(Nsc)
%  double precision MT_COEFF(Nsc),HT_COEFF,rhoc,cgmixc
%
%     do iu = 1,Nsc
%         nnu = (  (1/Molwt(iu))  + (1/Molwt(Nsc+1)) )**0.5
%         ddn = ( (nu(iu)**(0.333)) + (nu(Nsc+1)**(0.333)) )**2
%         Df(iu) = 0.1*(1.013d-02)*(Tchh**1.75)*(nnu/ddn)*(1/P)
%         MT_COEFF(iu) = Df(iu)*2.47/hyd_dia
%     enddo
%
% rhoc= 0.8                             !kg/m3
% cgmixc=1000                           !J/molK
% lambda_g = (2.66d-04)*(Tchh**0.805)
% HT_COEFF = 2.47*lambda_g/hyd_dia
%
% end subroutine MT_HT_COEF