
run simu_exac.m
close all
clear

load exact.mat
ord_num = 3;
ir_num  = 5;
n_RK    = 3;


figure_at_time=[0,0.1,0.3,0.5,1,20];

% UStore    = zeros((ord_num+1)*(10*2^(ir_num))+1,ir_num,ord_num);
% UexcStore = zeros((ord_num+1)*(10*2^(ir_num))+1,ir_num,ord_num);

L2_ErrorStore=zeros(ord_num,ir_num);
L2_OrderStore=zeros(ord_num,ir_num-1);
% LInf_ErrorStore=zeros(ord_num,ir_num);
% LInf_OrderStore=zeros(ord_num,ir_num-1);

fig = 0; 
for Ord = 1 : ord_num
    for ir = 1 : ir_num
        
        Time = 0;
        Nelm = 10*2^(ir-1);
        dx = period/Nelm;
        x = 0:dx:period;
        elm_size = Ord+1;
        
        dt = cfl * dx;
        Tsteps = floor((Tfinal-0.1*dt)/dt)+1;
        dt_final = Tfinal - (Tsteps-1) * dt;
        
        U0 = setInitial_smo(Nelm,elm_size,x);
        U = U0;
        plot_uh( U0,Ord,Nelm,x ,"exact")
        [ Amat,Pvmat,massMat,massMat_inv,mu_massMat ] = getAmat(Ord,Nelm,x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if Ord == ord_num && ir == ir_num
%             fig = fig+1;
%             figure(fig)
%             plot_uh( U,Ord,Nelm,x ,"exact")
%             grid on
%             xlabel('x')
%             ylabel('u')
%             Title_str = strcat('t=',num2str(Time));
%             title(Title_str);
%             % fig_name = 'smo_init';
%             % saveas(fig,['./smooth/eps/',fig_name],'eps')
%             % saveas(fig,['./smooth/fig/',fig_name],'fig')
%         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        for nt = 1:Tsteps
            if nt == Tsteps-1
                dt = dt_final;
            end
            U = RKn( Ord,x,Nelm,U,Amat,Pvmat,massMat,massMat_inv,mu_massMat,n_RK,dt,Time);
            Time = Time+dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%             if Ord == ord_num && ir == ir_num
%                 for p=1:size(figure_at_time,2)
%                     if figure_at_time(p)-dt/2 < Time && Time <= figure_at_time(p)+dt/2
%                         fig = fig+1;
%                         figure(fig)
%                         plot_uh( U,Ord,Nelm,x ,"numerical")
%                         hold on
%                         xlabel('x')
%                         ylabel('u')
%                         Title_str = strcat('t=',num2str(Time));
%                         title(Title_str);
%                         Tstr = strrep(num2str(Time),'.','_');
% %                         fig_name = strcat('sin_o',num2str(Ord),'i',num2str(ir),'t',Tstr);
% %                         saveas(fig,['./smooth/eps/',fig_name],'eps')
% %                         saveas(fig,['./smooth/fig/',fig_name],'fig')
%                     end
%                 end
%             end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        end
         Ue = get_exac( elm_size,Nelm,x,Ord_e,ir_e,Uexc,period );
%          L2_ErrorStore(Ord,ir)= L2err_integral( U,Ue,Ord,Nelm,x );
         L2_ErrorStore(Ord,ir)= l2err_discrete( U,Ue,Ord,Nelm);
    end
    L2_OrderStore(Ord,:)=L2_OrderStore(Ord,:)+ErrorOrder(L2_ErrorStore(Ord,:));
end
hold on
plot_uh( Ue,Ord,Nelm,x ,"exact")
plot_uh( U,Ord,Nelm,x ,"numerical")
    % fprintf('Infos:\nOrder_of_polynomial: %s \nMeshsize: %s \nRK_order: %s \nCFL: %s\n',...
    % num2str(ord_num),num2str(10*2^(ir_num-1)),num2str(n_RK),num2str(cfl))
format short
disp(L2_OrderStore)