clear; 
clc;

Tem = 5.8d2+2.73d2;
ini_temp = Tem;

time_beg = 0;
time_end = 2200;
delt = 1.0;
time_int = delt;


Nsc=5;
NstripsL2 = 9;
% NstripsL2 = 25;
NwgpL2= NstripsL2+1;
TNwgp =  NwgpL2;

NswL2_g = Nsc; % gas species  in bulk and in solid phase
NswL2_s = 1;   % surface species  in solid phase
NswL2= NswL2_g + NswL2_s; %total number of species

Twcs = Nsc + (NwgpL2*NswL2) + 2;
Nc= 10;
% Nc= 25;

yall=zeros(Twcs,Nc);

%% initial conditions:
%yall(Twcs-1,:)=ini_temp; %T_f
%yall(Twcs,:)=ini_temp;





%i_mol=[0.0, 0.1, 0.0, 0.1, 0.0, 0.8];    % at time 0sec
i_mol=[0.1, 0.2, 0.1, 0.2, 0.5, 0.7893];
%i_mol=[0.0, 0.1, 0.0, 0.1, 0.0, 0.8];

if (time_int >= 4.2 && time_int <= 14.1 )
    i_mol=[0.0, 0.1, 0.0, 0.1, 0.0107, 0.7893];  % 4.2 TO 14.1 sec
end

if (time_int >= 14.1 && time_int <= 24 )
    i_mol=[0.0, 0.1, 0.0, 0.1, 0.0, 0.8];        %14.1 to 24  sec
end

if (time_int >= 24 && time_int <= 34.2 )
    i_mol=[0.0045, 0.1, 0.0, 0.1, 0.0, 0.7955];  %24 to 34.3    sec
end

if (time_int >= 34.3 && time_int <= 38 )
    i_mol=[0.0, 0.1, 0.0, 0.1, 0.0, 0.8];        % 34.3 to 38    sec
end


for ii = 1:Nc
    k = 1;
    for i = 1:Nsc
        yall(k,ii)= i_mol(i);
        k = k+1;
    end
    for j = 1:NwgpL2
        for kk = 1:NswL2
            if (kk <= NswL2_g)
                yall(k,ii) = i_mol(kk);
            else
                yall(k,ii)= 1;
            end
            k = k+1;
        end
    end
    yall(k,ii) = Tem;
    k = k+1;
    yall(k,ii) = Tem;
    k = k+1;
end


%for t0=0:13    
    [t,y] = ode15s(@evaluate_rhs,[0,30],yall(:));
    plot_sol_3D_new(t,y,[1 130:10:320],1:5,Nc)
    
    %yall=y(end,:)';
    %yall = reshape(yall,Twcs,Nc);
    %yall(Twcs-1,1) = 5.8d2+2.73d2 + 0.05*(t0+1);
%end

% while (time_int <= time_end)
%     yall(Twcs-1,1) = ini_temp  + (0.05*time_beg); % Set Dirichlet value
%     %CTotal = P/(Rconst*yall(Twcs-1,1));
%     
%     %Implicit Euler:
%     yall = yall + delt*evaluate_rhs(yall);
%     
%     yall=abs(yall);
%     
%     time_beg = time_int;
%     time_int = time_int + delt;
% end