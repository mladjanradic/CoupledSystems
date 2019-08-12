% yall = alle U zum Zeitpunkt k
% Berechnet werden soll Zeitpunkt k+1

function [yall_new] = evaluate_rhs(t,yall)

t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants we need
V = 1.65e0;
tort = 3.0;
dp   = 1e-7;
Molwt=[31.99886, 18.01528, 2.015894, 44.00964, 28.01021, 28.01344];
nu=[ 16.6 , 12.7 , 7.07 , 26.9,  18.9,  17.9 ];
x0=0;
xf=7.62e-2;
%dh = 1.2e-3;
sub_thick = 150.0e-06;
thick_L2 = 7.2e-5;
Dch =  1.e-3; %(2*sub_thick) + (2*thick_L2) + dh
dh = Dch - ((2*sub_thick) + (2*thick_L2) );
CPSI = 600; %(0.0254/Dch)**2;
Ageo = 4*CPSI*dh/(0.0254^2);     % geometric surface area
epsilong = CPSI*dh*dh/(0.0254^2);  % %   volume fraction of monolith
thick_w = (sub_thick) +thick_L2;
hyd_dia = dh;
Area = 3.14*dh*dh/4;
Volume = Area*xf;
SV = 51430;   % hr-1
% V = SV*xf/3600
gama = 76.7; % storage capacity in mol/m3
Nsc= 5;      % channel species number

Nc= 10;      % channel grid points
NstripsL2 = 9;

%Nc= 25;      % channel grid points
%NstripsL2 = 25;

NwgpL1 = 0;
NwgpL2= NstripsL2+1;
TNwgp =  NwgpL2; % total no of gp in both layers
NswL2_g = Nsc; % gas species  in bulk and in solid phase
NswL2_s = 1;   % surface species  in solid phase
NswL2= NswL2_g + NswL2_s; %total number of species
num_spe = NswL2;
Twcs = Nsc+1+NswL2+1; % no of eqns in one cascade
no_rxns =6;
thick_L1 = 0.0;
delz = (xf-x0)/(Nc-1);
delwc_L2 = thick_L2/NstripsL2;
adtemp = 5.8d2+2.73d2;
Tem =  5.8d2+2.73d2 ;
ini_temp = Tem;
solid_ini_temp = Tem;
Rconst = 8.314d0;
P = 1.0d5;
%epsilong =  0.7d0;
epsilonwc = 0.41d0;
%epsilonSolid = 0.4d0
% substrate
den_solid = 1800; % kg/m3
% Washcoat
den_WC = 1000;
cp_WC = 1000;
Twcs = Nsc + (NwgpL2*NswL2) + 2;
WCHEAT_COND=1.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Reshape it, because matlab can't handle the "matrix-form"
yall = reshape(yall,Twcs,Nc);

%Set Dirichlet value
%yall(Twcs-1,1) = 5.8d2+2.73d2 + 0.05*t;

for chan_pos=1:Nc


U = yall(:,chan_pos);
term1=zeros(size(U));
term2=zeros(size(U));
term3=zeros(size(U));
DELTA=zeros(size(U));

count = 1;

Tchh = U(Twcs-1);
[MT_COEF,HT_COEF,rho,cgmix] = MT_HT_COEF(Tchh);

%% channel species
for i = 1:Nsc
    if (chan_pos == 1)
        term1(count) = 0;
        term2(count) = 0;
        term3(count) = 0;
        DELTA(count) = 0;  % here the 1st gridf point in the bulk phase is assumed to be the inlet point
        count = count+1;   %that is the reason the derivative is assumed to be zero as there won't be any change in the concentations
    else
        term1(count) = -(V/delz)*(U(count)- yall(count,chan_pos-1));
        term2(count) = -MT_COEF(i)*Ageo*( U(count)-U(count+Nsc) )/epsilong;
        term3(count) = 0;
        DELTA(count) = term1(count)+term2(count)+term3(count);
        count = count+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% washcoat species %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_spe = NswL2;
num_gas_spe = NswL2_g;
num_surf_spe = NswL2_s;
delwc_L = delwc_L2;
no_rxn = NswL2_g+NswL2_s;

for j = 1:NwgpL2
    i_mf=zeros(num_spe,1);
    for i =1:num_spe
        i_mf(i) = U(count+i-1);
    end
    for i=1:Nsc
        Dknud(i) = dp*97*( sqrt(U(Twcs)/Molwt(i)) );
    end
    Tchh = U(Twcs);
    
    for i=1:Nsc
        nnu = (  (1/Molwt(i))  + (1/Molwt(Nsc+1)) )^0.5;
        ddn = ( (nu(i)^(0.333)) + (nu(Nsc+1)^(0.333)) )^2;
        Dext(i)=(1.013d-02)*(Tchh^1.75)*(nnu/ddn)*(1/P);
        Deff(i) = (epsilonwc/tort)*Dext(i)*Dknud(i)/(Dext(i)+Dknud(i));
    end
    [H_Rxn,r_fm,i_rs] = rate_calc_one(i_mf,Tchh);
    for i = 1:no_rxn
        Heat_Rxn(j,i) = H_Rxn(i);
        indivi_rates(j,i) = i_rs(i);
    end
    
    % Boundary or inner points:
    if (j==1 )
        condition = 1; % at the gas/solid interface
    elseif (j == NwgpL2)
        condition = 3;
    else
        condition = 2;
    end
    
    for kk = 1:NswL2
        if (kk > num_gas_spe  && kk <= num_spe )
            term1(count) = r_fm(kk);
            term2(count) = 0;
            term3(count) = 0;
            DELTA(count) = term1(count)+term2(count)+term3(count);
            count = count +1;
        else
            if (condition == 1) % at the gas/solid interface
                cflux1 = 2*(U(count+NswL2)-U(count))/delwc_L;
                cflux2 = 2*MT_COEF(kk)*( U(count)-U(kk))/Deff(kk);
                factor1 = (Deff(kk))/delwc_L;
            elseif(condition == 2) % internal points
                cflux1 = (U(count+NswL2)-U(count))/delwc_L;
                cflux2 = (U(count)-U(count-NswL2))/delwc_L;
                factor1 = (Deff(kk))/(delwc_L);
            else  %(condition == 3) % for points at the washcoat end
                cflux1 = 2*(U(count)-U(count))/delwc_L;
                cflux2 = 2*(U(count)-U(count-NswL2))/delwc_L;
                factor1 = (Deff(kk))/(delwc_L);
            end
            CTotal = P/(Rconst*U(Twcs));
            term1(count) = factor1 *(cflux1-cflux2)/epsilonwc;
            term2(count) = (gama*r_fm(kk)/(CTotal*epsilonwc));
            term3(count) = 0;
            DELTA(count) = term1(count)+term2(count)+term3(count);
            count = count+1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (chan_pos == 1)
    U((chan_pos*Twcs)-1) = ini_temp;
    term1(count) = 0;
    term2(count) = 0;
    term3(count) = 0;
    DELTA(count) = 0;
    count = count+1;
else
    term1(count) = -(V/(delz))*(U(Twcs-1)-yall(Twcs-1,chan_pos-1));
    term2(count) = (HT_COEF*Ageo*(U(Twcs) - U(Twcs-1))/(rho*cgmix*epsilong));
    term3(count) = 0;
    DELTA(count) = term1(count)+term2(count)+term3(count);
    count = count+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% solid energy balance%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(chan_pos == 1)
    condi = 1;
elseif (chan_pos == Nc)
    condi = 3;
else
    condi = 2;
end


denom1 = (den_WC*cp_WC*(1-epsilong));
denom2 = 0;
denominator = denom1 + denom2;
tfactor1 = (WCHEAT_COND*(1-epsilong));


if (condi == 1) % at the gas/solid interface
    tflux1 = 2*( yall(Twcs,(chan_pos+1))-U(Twcs) )/delz;
    tflux2 = 2*( U(Twcs)-U(Twcs) )/delz;
elseif (condi == 2) % internal points
    tflux1 = ( yall(Twcs,(chan_pos+1)) -U(Twcs) )/delz;
    tflux2 = ( U(Twcs) -yall(Twcs,(chan_pos-1)) )/delz;
else % (condition == 3) then  %for all internal points
    tflux1 = 2*( U(Twcs) -U(Twcs) )/delz;    % at the solid/solid interface
    tflux2 = 2*( U(Twcs) -yall(Twcs,(chan_pos-1)) )/delz;
end
term1(count) = (tfactor1/(denominator*delz))*(tflux1 - tflux2);
term3(count) = -(HT_COEF*Ageo*( U(Twcs) - U(Twcs-1)))/denominator;
sumTot= 0;
for i = 1:no_rxn  % here it for number of gas species, neglecting N2
    sum = 0;
    sum1=zeros(TNwgp,1);
    
    for ig = 1:TNwgp
        sum1(ig) = (indivi_rates(ig,i)*Heat_Rxn(ig,i));
    end
    
    sum = ( sum1(1)+sum1(TNwgp) );
    
    for ig = 2:TNwgp-1;
        if (mod(ig,2)==0)
            sum = ( 4*sum1(ig) ) + sum;
        else
            sum = ( 2*sum1(ig) ) + sum;
        end
    end
    sumTot = (sum*thick_L2/(NstripsL2*3) ) + sumTot;
end

term2(count) = -(gama*sumTot*Ageo/(denominator));
DELTA(count) = term1(count)+term2(count)+term3(count);


yall_new(:,chan_pos) = DELTA;
end %%For-Schleife über z-Achse

yall_new=yall_new(:);

