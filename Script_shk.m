clear
close all

Tfinal =0;

figure_at_time=[0,1,2,5,Tfinal];
%  figure_at_time=[0,1,3,5,10,20];
ord_num = 1;
ir_num = 5;
cfl = 0.1;
n_RK  = 3;
period = 1;
% indicator of the initial data
% 1:one shock; 2:two shocks; 3:three shocks
CS = 3;

P0 = zeros(CS,1);
Q0 = zeros(CS,1);
S0 = zeros(CS,1);
switch CS
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
        P0 = [1; 0.8; 0.12];
        Q0 = [0.1; 0.5; 0.8];
        S0 = [0.7; 0.4; 0.2];
end

UStore    = zeros((ord_num+1)*(10*2^(ir_num))+1,ir_num,ord_num);
UexcStore = zeros((ord_num+1)*(10*2^(ir_num))+1,ir_num,ord_num);

fig = 0;
for Ord = ord_num:ord_num
    for ir = ir_num:ir_num
        
        Time = 0;        
        Nelm = 10*2^(ir-1);
        dx = period/Nelm;
        x = 0:dx:period;
        elm_size = Ord+1;
        
        dt = cfl * dx;
        Tsteps = floor((Tfinal-0.1*dt)/dt)+1;
        dt_final = Tfinal - (Tsteps-1) * dt;
        
        TV = zeros(Tsteps);
        U0 = setInitial_shock(Nelm,elm_size,x,CS,period,P0,Q0,S0);
        
        [ Amat,Pvmat,massMat,massMat_inv,mu_massMat ] = getAmat(Ord,Nelm,x);
        U = U0;
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fig = fig+1;
        figure(fig)
        plot_uh( U,Ord,Nelm,x ,"exact")
        grid on
        xlabel('x')
        ylabel('u')
        Title_str = strcat('t=',num2str(Time));
        title(Title_str);
        fig_name = strcat('sk',num2str(CS),'_init');
        saveas(fig,['./simulations/eps/',fig_name],'eps')
        saveas(fig,['./simulations/fig/',fig_name],'fig')
        
        for nt = 1:Tsteps
            
            if nt == Tsteps-1
                dt = dt_final;
            end
            
            if Ord == 0
                U = RKn( Ord,x,Nelm,U,Amat,Pvmat,massMat,massMat_inv,mu_massMat,n_RK,dt,Time);
            else
                U = RKn_limiter( Ord,x,Nelm,U,Amat,Pvmat,massMat,massMat_inv,mu_massMat,n_RK,dt,Time,P0 );
            end
            
%             TV(nt) = total_variation(U,Nelm,Ord,x);
            Time = Time+dt;
            
            for p=1:size(figure_at_time,2)
                if figure_at_time(p)-dt/2 < Time && Time <= figure_at_time(p)+dt/2
                    fig = fig+1;
                    figure(fig)
                    plot_uh( U,Ord,Nelm,x ,"numerical")
                    hold on
                    [ Uexc,~,~ ] = shock_pkns_solu( P0,Q0,S0,Nelm,elm_size, x ,period,Time,CS );
                    plot_uh( Uexc,Ord,Nelm,x ,"exact")
                    grid on
                    xlabel('x')
                    ylabel('u')
                    Title_str = strcat('t=',num2str(Time));
                    title(Title_str);
                    %                    legend('LDG','Exact')
                    Tstr = strrep(num2str(Time),".","_");
                    fig_name = strcat('sk',num2str(CS),'_o',num2str(Ord),'i',num2str(ir),'t',Tstr);
                    saveas(fig,['./simulations/eps/',fig_name],'eps')
                    saveas(fig,['./simulations/fig/',fig_name],'fig')
                    
                end
            end
        end
    end
end
fprintf('Infos:\nOrder_of_polynomial: %s \nMeshsize: %s \nshocks: %s \nRK_order: %s \nCFL: %s\n',...
    num2str(ord_num),num2str(10*2^(ir_num-1)),num2str(CS),num2str(n_RK),num2str(cfl))






