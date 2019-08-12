function [] = plot_sol_3D_new(t,y,T,species,chan_pos)

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

Z=linspace(x0,xf,Nc);
R=linspace(0,thick_L2,NwgpL2);

[XTZ,YTZ] = meshgrid(t(T),Z);
[XTR,YTR] = meshgrid(t(T),R);
Xf = ['O2 ','H20','H2 ','CO2','CO '];


for j=species
    if j<6
        counter=1;
        ynew=y(T,j:Twcs:Twcs*Nc);
        figure(j)
        subplot(1,1+length(chan_pos),counter)
        plot3(XTZ,YTZ,ynew,'b','LineWidth',2);
        set(gca,'FontSize',14); 
        TITLE=['Specie ' num2str(j) ' of X_f - ' Xf(3*(j-1)+1:3*j)];
        title(TITLE);
        xlabel('Time'), ylabel('Channel'), zlabel(Xf(3*(j-1)+1:3*j))
        
        for k=chan_pos
            counter= counter+1;
            mystart= Twcs*(k-1)+Nsc+j;
            mystop = Twcs*k-2;
            ynew=y(T,mystart:6:mystop);
            subplot(1,1+length(chan_pos),counter)
            plot3(XTR,YTR,ynew,'b','LineWidth',2)
            set(gca,'FontSize',14); 
            TITLE=['Specie ' num2str(j) ' of X_s - ' Xf(3*(j-1)+1:3*j)];
            title(TITLE);
            %xlabel('Time'), ylabel(['Washcoat Channel Position: ' num2str(k)]), zlabel(Xf(3*(j-1)+1:3*j))
            xlabel('Time'), ylabel(['Washcoat Outlet']), zlabel(Xf(3*(j-1)+1:3*j))
        end 
    end
    
    if j==6
        figure(j)
        counter=1;
        for k=chan_pos
            mystart= Twcs*(k-1)+Nsc+j;
            mystop = Twcs*k-2;
            ynew=y(T,mystart:6:mystop);
            subplot(1,length(chan_pos),counter)
            plot3(XTR,YTR,ynew,'b','LineWidth',2)
            set(gca,'FontSize',14);
            TITLE=['Specie ' num2str(j) ' - \theta' ];
            title(TITLE);
            xlabel('Time'), ylabel('Washcoat Outlet'), zlabel('\theta')
            counter=counter+1;
        end   
    end
    
    if j==7
        ynew=y(T,Twcs-1:Twcs:Twcs*Nc);
        figure(j)
        subplot(1,2,1)
        plot3(XTZ,YTZ,ynew,'b','LineWidth',2)
        set(gca,'FontSize',14);
        TITLE='Bulk Energy T_f';
        title(TITLE);
        xlabel('Time'), ylabel('Channel'), zlabel('T_f')
        
        ynew=y(T,Twcs:Twcs:Twcs*Nc);
        subplot(1,2,2)
        plot3(XTZ,YTZ,ynew,'b','LineWidth',2)
        set(gca,'FontSize',14);
        TITLE='Solid Energy T_s';
        title(TITLE);
        xlabel('Time'), ylabel('Channel'), zlabel('T_s')
        
    end
    
end

end

