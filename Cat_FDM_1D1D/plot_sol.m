function [] = plot_sol(t,y,T,species)

%[T L] = size(y);

x0=0;
xf=7.62e-2;
Nsc= 5;      % channel species number
Nc= 10;      % channel grid points
NstripsL2 = 9;
delz = (xf-x0)/(Nc-1);
thick_L2 = 7.2e-5;
delwc_L2 = thick_L2/NstripsL2;
NwgpL2= NstripsL2+1;
NswL2_g = Nsc;
NswL2_s = 1;
NswL2= NswL2_g + NswL2_s;
Twcs = Nsc + (NwgpL2*NswL2) + 2;

R=linspace(0,7.2e-5,NstripsL2+1);
Z=linspace(x0,xf,Nc);
[X Y] = meshgrid(Z,R);

Xf = ['O2 ','H20','H2 ','CO2','CO '];

for i=T
    yT = y(i,:);
    yT=reshape(yT,Twcs,Nc);
    
    for j=species
        
        if(j<6)
            figure(j)
            subplot(1,2,1)
            
            plot(x0:delz:xf,yT(j,:))
            TITLE=['Specie ' num2str(j) ' of X_f - ' Xf(3*(j-1)+1:3*j) ' - Timestep ' num2str(t(i))];
            title(TITLE);
            xlabel('Channel'), ylabel(['Specie ' num2str(j)])
            
            subplot(1,2,2)
            plot3(X,Y,yT(5+j:6:Twcs-2,:)')
            TITLE=['Specie ' num2str(j) ' of X_s - ' Xf(3*(j-1)+1:3*j) ' - Timestep ' num2str(t(i))];
            title(TITLE);
            ylabel('Washcoat'), xlabel('Channel'), zlabel(['Specie ' num2str(j)])
        end
        
        if(j==6)
            figure(j)
            plot3(X,Y,yT(8:6:Twcs-2,:))
            TITLE=['6th Specie of X_s - Theta - Timestep ' num2str(t(i))];
            title(TITLE);
            xlabel('Washcoat'), ylabel('Theta')
        end
        
        
        if (j==7)
            figure(j)
            subplot(1,2,1)
            plot(x0:delz:xf,yT(Twcs-1,:))
            TITLE=['Bulk Energy T_f - Timestep ' num2str(t(i))];
            title(TITLE);
            xlabel('Channel'), ylabel('T_f')

            subplot(1,2,2)
            plot(x0:delz:xf,yT(Twcs,:))
            TITLE=['Solid Energy T_s - Timestep ' num2str(t(i))];
            title(TITLE);
            xlabel('Channel'), ylabel('T_f')
        end
    end
    
    if(length(T)>1)
        pause
    end
    
end





% function [] = plot_sol(t,y)
%
% [T L] = size(y);
%
% x0=0;
% xf=7.62e-2;
% Nsc= 5;      % channel species number
% Nc= 10;      % channel grid points
% NstripsL2 = 9;
% delz = (xf-x0)/(Nc-1);
% thick_L2 = 7.2e-5;
% delwc_L2 = thick_L2/NstripsL2;
% NwgpL2= NstripsL2+1;
% NswL2_g = Nsc;
% NswL2_s = 1;
% NswL2= NswL2_g + NswL2_s;
% Twcs = Nsc + (NwgpL2*NswL2) + 2;
%
% for i=T:T
%     yT = y(i,:);
%     yT=reshape(yT,Twcs,Nc);
%
% figure(1)
%
% subplot(5,1,1)
% plot(x0:delz:xf,yT(1,:))
% TITLE=['1st Specie of X_f - O2 - Timestep ' num2str(t(i))];
% title(TITLE);
% xlabel('Channel'), ylabel('O2')
%
% subplot(5,1,2)
% plot(x0:delz:xf,yT(2,:))
% TITLE=['2nd Specie of X_f - H2O - Timestep ' num2str(t(i))];
% title(TITLE);
% xlabel('Channel'), ylabel('H2O')
%
% subplot(5,1,3)
% plot(x0:delz:xf,yT(3,:))
% TITLE=['3rd Specie of X_f - H2 - Timestep ' num2str(t(i))];
% title(TITLE);
% xlabel('Channel'), ylabel('H2')
%
% subplot(5,1,4)
% plot(x0:delz:xf,yT(4,:))
% TITLE=['4th Specie of X_f - CO2 - Timestep ' num2str(t(i))];
% title(TITLE);
% xlabel('Channel'), ylabel('CO2')
%
% subplot(5,1,5)
% plot(x0:delz:xf,yT(5,:))
% TITLE=['5th Specie of X_f - CO - Timestep ' num2str(t(i))];
% title(TITLE);
% xlabel('Channel'), ylabel('CO')
%
%
% R=linspace(0,7.2e-5,NstripsL2+1);
% Z=x0:delz:xf;
%
% % for j=1:Nc
% %     figure(j+1)
% %
% %     subplot(6,1,1)
% %     plot(R,yT(6:6:Twcs-2,j))
% %     TITLE=['1st Specie of X_s - O2 - Timestep ' num2str(t(i)) ' Channel position: z=' num2str(Z(j))];
% %     title(TITLE);
% %     xlabel('Washcoat'), ylabel('O2')
% %
% %     subplot(6,1,2)
% %     plot(R,yT(7:6:Twcs-2,j))
% %     TITLE=['2nd Specie of X_s - H20 - Timestep ' num2str(t(i)) ' Channel position: z=' num2str(Z(j))];
% %     title(TITLE);
% %     xlabel('Washcoat'), ylabel('H20')
% %
% %     subplot(6,1,3)
% %     plot(R,yT(8:6:Twcs-2,j))
% %     TITLE=['3rd Specie of X_s - H2 - Timestep ' num2str(t(i)) ' Channel position: z=' num2str(Z(j))];
% %     title(TITLE);
% %     xlabel('Washcoat'), ylabel('H2')
% %
% %     subplot(6,1,4)
% %     plot(R,yT(8:6:Twcs-2,j))
% %     TITLE=['4th Specie of X_s - CO - Timestep ' num2str(t(i)) ' Channel position: z=' num2str(Z(j))];
% %     title(TITLE);
% %     xlabel('Washcoat'), ylabel('CO')
% %
% %     subplot(6,1,5)
% %     plot(R,yT(8:6:Twcs-2,j))
% %     TITLE=['5th Specie of X_s - CO2 - Timestep ' num2str(t(i)) ' Channel position: z=' num2str(Z(j))];
% %     title(TITLE);
% %     xlabel('Washcoat'), ylabel('CO2')
% %
% %     subplot(6,1,6)
% %     plot(R,yT(8:6:Twcs-2,j))
% %     TITLE=['6th Specie of X_s - Theta - Timestep ' num2str(t(i)) ' Channel position: z=' num2str(Z(j))];
% %     title(TITLE);
% %     xlabel('Washcoat'), ylabel('Theta')
% %
% % end
%
%
%     [X Y] = meshgrid(R,Z);
%
%     figure(2)
%
%     subplot(6,1,1)
%     plot3(X,Y,yT(6:6:Twcs-2,:))
%     TITLE=['1st Specie of X_s - O2 - Timestep ' num2str(t(i))];
%     title(TITLE);
%     xlabel('Washcoat'), ylabel('O2')
%
%     subplot(6,1,2)
%     plot3(X,Y,yT(7:6:Twcs-2,:))
%     TITLE=['2nd Specie of X_s - H20 - Timestep ' num2str(t(i))];
%     title(TITLE);
%     xlabel('Washcoat'), ylabel('H20')
%
%     subplot(6,1,3)
%     plot3(X,Y,yT(8:6:Twcs-2,:))
%     TITLE=['3rd Specie of X_s - H2 - Timestep ' num2str(t(i))];
%     title(TITLE);
%     xlabel('Washcoat'), ylabel('H2')
%
%     subplot(6,1,4)
%     plot3(X,Y,yT(8:6:Twcs-2,:))
%     TITLE=['4th Specie of X_s - CO - Timestep ' num2str(t(i))];
%     title(TITLE);
%     xlabel('Washcoat'), ylabel('CO')
%
%     subplot(6,1,5)
%     plot3(X,Y,yT(8:6:Twcs-2,:))
%     TITLE=['5th Specie of X_s - CO2 - Timestep ' num2str(t(i))];
%     title(TITLE);
%     xlabel('Washcoat'), ylabel('CO2')
%
%     subplot(6,1,6)
%     plot3(X,Y,yT(8:6:Twcs-2,:))
%     TITLE=['6th Specie of X_s - Theta - Timestep ' num2str(t(i))];
%     title(TITLE);
%     xlabel('Washcoat'), ylabel('Theta')
%
%
%
% figure(Nc+2)
% subplot(2,1,1)
% plot(x0:delz:xf,yT(Twcs-1,:))
%
% TITLE=['Bulk Energy T_f - Timestep ' num2str(t(i))];
% title(TITLE);
% xlabel('Channel'), ylabel('T_f')
%
% subplot(2,1,2)
% plot(x0:delz:xf,yT(Twcs,:))
% TITLE=['Solid Energy T_s - Timestep ' num2str(t(i))];
% title(TITLE);
% xlabel('Channel'), ylabel('T_f')
% end
% end