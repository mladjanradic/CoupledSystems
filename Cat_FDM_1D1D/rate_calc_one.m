function [H_Rxn,r_f,rate] = rate_calc_one(i_m,Temp)

%% Everything we need for this function:
Nsc=5;
Rconst = 8.314d0;
P = 1.0d5;
NswL2_g = Nsc; % gas species  in bulk and in solid phase
NswL2_s = 1;   % surface species  in solid phase
NswL2= NswL2_g + NswL2_s; %total number of species 


Cg=zeros(NswL2,1);
H_Rxn = zeros(NswL2,1);
rate=zeros(NswL2,1);
r_f=zeros(NswL2,1);


for igt = 1:Nsc+1
    Cg(igt)=i_m(igt);
end

Ce2O3_O = Cg(6);
Ce2O3_ = 1-Ce2O3_O;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cg_O2  = Cg(1)*P/(Rconst*Temp);
Cg_H2O = Cg(2)*P/(Rconst*Temp);
Cg_H2  = Cg(3)*P/(Rconst*Temp);
Cg_CO2 = Cg(4)*P/(Rconst*Temp);
Cg_CO  = Cg(5)*P/(Rconst*Temp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumCg = P/(Rconst*Temp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_CO    = -121287.3305583559;
H_H2    = -9226.295742641230;
H_O2    = -11222.99746307711;
H_CO2   = -412764.0956540055;
H_H2O   = -255701.4820603337;
S_CO    =  183.9187223864227;
S_H2    =  118.0583120890414;
S_O2    =  189.8424309558777;
S_CO2   =  189.2775313347509;
S_H2O   =  172.7032053657214;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_Ce2O4c        = -2.896594e+005;
fH_Ce2O4Ce2O4   = 5.434121e+004;
S_Ce2O4         = 5.029200e+000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_Ce2O4 = H_Ce2O4c + fH_Ce2O4Ce2O4 * (Ce2O3_O)^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A1      =  1.605787e+004;
E1      =  1.000000e+004;
A2      =  1.993722e+014;
E2      =  1.607071e+005;
A3      =  4.099917e+014;
E3      =  2.301590e+005;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A1b     = A1 * exp((-2 * S_Ce2O4 + S_O2)/Rconst) * sumCg;
E1b     = E1 - 2 * H_Ce2O4 + H_O2;
A2b     = A2 * exp(-(S_Ce2O4 + S_H2 - S_H2O)/Rconst);
E2b     = E2 - H_Ce2O4 - H_H2 + H_H2O;
A3b     = A3 * exp(-(S_Ce2O4 + S_CO - S_CO2)/Rconst);
E3b     = E3 - H_Ce2O4 - H_CO + H_CO2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Oxidation of Cer by O2
rate(1) = A1  * exp(-E1 /(Rconst*Temp)) * Cg_O2 * (Ce2O3_)^2;
rate(2) = A1b * exp(-E1b/(Rconst*Temp)) * (Ce2O3_O)^2;
% Oxidation of Cer by H2O
rate(3) = A2  * exp(-E2 /(Rconst*Temp)) * Cg_H2O * (Ce2O3_);
rate(4) = A2b * exp(-E2b/(Rconst*Temp)) * Cg_H2  * (Ce2O3_O);
% Oxidation of Cer by CO2
rate(5) = A3  * exp(-E3 /(Rconst*Temp)) * Cg_CO2 * (Ce2O3_);
rate(6) = A3b * exp(-E3b/(Rconst*Temp)) * Cg_CO  * (Ce2O3_O);
r_f(1) = -rate(1) + rate(2);
r_f(2) = -rate(3) + rate(4);
r_f(3) =  rate(3) - rate(4);
r_f(4) = -rate(5) + rate(6);
r_f(5) =  rate(5) - rate(6);
r_f(6) =  ((2*rate(1))-(2*rate(2))+rate(3) ...
            -rate(4)+rate(5)-rate(6));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_CO    = -121287.3305583559;
H_H2    = -9226.295742641230;
H_O2    = -11222.99746307711;
H_CO2   = -412764.0956540055;
H_H2O   = -255701.4820603337;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_Ce2O4c        = -2.896594e+005;
fH_Ce2O4Ce2O4   = 5.434121e+004;
H_Ce2O4 = H_Ce2O4c + fH_Ce2O4Ce2O4 * (Ce2O3_O)^2;
H_Ce2O3 = 0;

H_Rxn(1) = 2* H_Ce2O4-(H_O2+(2*H_Ce2O3));
H_Rxn(2) = -H_Rxn(1);
H_Rxn(3) = H_H2 + H_Ce2O4-(H_H2O +H_Ce2O3);
H_Rxn(4) = -H_Rxn(3);
H_Rxn(5) = H_CO + H_Ce2O4-(H_CO2 +H_Ce2O3);
H_Rxn(6) = -H_Rxn(5);


end