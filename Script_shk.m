clear
close all

Tfinal = 1;

figure_at_time=[0,0.5,1,3,5];
fig_Ey_at_time = 0:0.1:(Tfinal);

ord_num = 2;
ir_num  = 4;
cfl     = 0.1;
period  = 1;
% RK3 is better for shock solution
n_RK    = 3;
shock_type = 1;


P0 = zeros(shock_type,1);
Q0 = zeros(shock_type,1);
S0 = zeros(shock_type,1);
switch shock_type
    % Q-location of shocks; P-speed of transport; S-size of shock at Q
    case 1
        P0 = 0.333;
        Q0 = 0.1;
        S0 = 0.1;
    case 2
        P0 = [0.3; 0.1];
        Q0 = [0.2; 0.5];
        S0 = [0.4; 0.2];
    case 3
        P0 = [1  ; 0.8; 0.12];
        Q0 = [0.1; 0.5; 0.8];
        S0 = [0.7; 0.4; 0.2];
end

UStore    = zeros((ord_num+1)*(10*2^(ir_num))+1,ir_num,ord_num);
UexcStore = zeros((ord_num+1)*(10*2^(ir_num))+1,ir_num,ord_num);

Ey  = zeros(size(fig_Ey_at_time,2),2);
Eyc = zeros(size(fig_Ey_at_time,2),2);
Eyd = zeros(size(fig_Ey_at_time,2),2);

fig = 0;
nEy = 0;
for Ord = ord_num:ord_num
    for ir = ir_num:ir_num
        
        elm_size = Ord+1;
        % number of cells must be odd!!!
        Nelm = 10*2^(ir-1)+1;
        dx   = period/Nelm;
        x    = 0:dx:period;
            
        Time = 0;
        dt = cfl * dx;
        Tsteps = floor((Tfinal-0.1*dt)/dt)+1;
        dt_final = Tfinal - (Tsteps-1) * dt;
        %------------------------------------------------------------------
        U0 = setInitial_shock(Nelm,elm_size,x,shock_type,period,P0,Q0,S0);
        Uc  = U0;
        Ud  = U0;
        %------------------------------------------------------------------
        nEy = nEy+1;
        Eyc(nEy,:) = getEnergy( U0,Ord,Nelm,x);
        Eyd(nEy,:) = Eyc(nEy,:);
        %------------------------------------------------------------------
        fig = fig + 1;
        figure(fig)
        plot_uh( U0,Ord,Nelm,x ,"exact")
        grid on
        xlabel('x')
        ylabel('u')
        Title_str = strcat('t=',num2str(Time));
        title(Title_str);
        %------------------------------------------------------------------
        [ Amat_c,Pvmat_c,Pqmat_c,massMat_c,massMat_inv_c,mu_massMat_c ]...
            = getAmat(Ord,Nelm,x,'C','C');
        [ Amat_d,Pvmat_d,Pqmat_d,massMat_d,massMat_inv_d,mu_massMat_d ]...
            = getAmat(Ord,Nelm,x,'R','L');
        %------------------------------------------------------------------
        for nt = 1:Tsteps
            
            if nt == Tsteps
                dt = dt_final;
            end
            %--------------------------------------------------------------
            if Ord == 0
                Uc = RKn( Ord,x,Nelm,Uc,Amat_c,Pvmat_c,Pqmat_c,massMat_c,...
                    massMat_inv_c,mu_massMat_c,n_RK,dt,Time,'Csv' );
                Ud = RKn( Ord,x,Nelm,Ud,Amat_d,Pvmat_d,Pqmat_d,massMat_d,...
                    massMat_inv_d,mu_massMat_d,n_RK,dt,Time,'Dsp' ); 
            else
                Uc = RKn_limiter( Ord,x,Nelm,Uc,Amat_c,Pvmat_c,Pqmat_c,...
                    massMat_c,massMat_inv_c,mu_massMat_c,n_RK,dt,Time,'Csv');
                Ud = RKn_limiter( Ord,x,Nelm,Ud,Amat_d,Pvmat_d,Pqmat_d,...
                    massMat_d,massMat_inv_d,mu_massMat_d,n_RK,dt,Time,'Dsp');
            end
            %--------------------------------------------------------------            
            Time = Time+dt;
            %--------------------------------------------------------------
            for p=1:size(fig_Ey_at_time,2)
                if fig_Ey_at_time(p)-dt/2 < Time && Time <= fig_Ey_at_time(p)+dt/2
                    nEy = nEy+1;
                    Eyc(nEy,:) = getEnergy( Uc,Ord,Nelm,x);
                    Eyd(nEy,:) = getEnergy( Ud,Ord,Nelm,x);
                end
            end
            %--------------------------------------------------------------
            for p=1:size(figure_at_time,2)
                if figure_at_time(p)-dt/2 < Time && Time <= figure_at_time(p)+dt/2
                    fig = fig+1;
                    figure(fig)
                    plot_uh( Uc,Ord,Nelm,x,"numerical")
                    hold on
                    plot_uh( Ud,Ord,Nelm,x,"other")
                    [ Uexc,~,~ ] = shock_pkns_solu( P0,Q0,S0,Nelm,elm_size,...
                        x,period,Time,shock_type );
                    plot_uh( Uexc,Ord,Nelm,x ,"exact")
                    grid on
                    xlabel('x')
                    ylabel('u')
                    Title_str = strcat('t=',num2str(Time));
                    title(Title_str);
                    legend('LDG(csv)','LDG(dsp)','Exact')
                    
                end
            end
        end
    end
end

len = size(fig_Ey_at_time,2);
diff_Eyc0 = abs(Eyc(2:len,1)'-Eyc(1,1)*ones(1,len-1));
diff_Eyd0 = abs(Eyd(2:len,1)'-Eyd(1,1)*ones(1,len-1));
diff_Eyc1 = abs(Eyc(2:len,2)'-Eyc(1,2)*ones(1,len-1));
diff_Eyd1 = abs(Eyd(2:len,2)'-Eyd(1,2)*ones(1,len-1));

fig = fig+1;
figure(fig)
semilogy(fig_Ey_at_time(2:len),diff_Eyc0,'-^','LineWidth',1.5)
hold on
semilogy(fig_Ey_at_time(2:len),diff_Eyd0,'-x','LineWidth',1.5)
semilogy(fig_Ey_at_time(2:len),diff_Eyc1,'-v','LineWidth',1.5)
semilogy(fig_Ey_at_time(2:len),diff_Eyd1,'-o','LineWidth',1.5)
grid on

xlabel('t')
legend('E_0(csv)','E_0(dsp)','E_1(csv)','E_1(dsp)', 'Location','NorthEastOutside')
title('|E(t)-E(0)|')
