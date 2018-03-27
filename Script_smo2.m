% smooth travelling wave with initial condition from gCH equation.

close all
clear
load exact2.mat

ord_num = 3;
ir_num  = 5;
n_RK    = 4; 

L2_ErrorStore   = zeros(ord_num+1,ir_num);
L2_OrderStore   = zeros(ord_num+1,ir_num-1);
Linf_ErrorStore = zeros(ord_num+1,ir_num);
Linf_OrderStore = zeros(ord_num+1,ir_num-1);

fig = 0; 
for Ord = 0 : ord_num
    for ir = 1 : ir_num
        
        % Nelm must be odd.
        Nelm = 10*2^(ir-1)+1;
        dx = period/Nelm;
        x = -period/2:dx:period/2;
        elm_size = Ord+1;
        
        Time = 0;        
        dt = cfl * dx;
        Tsteps = floor((Tfinal-0.1*dt)/dt)+1;
        dt_final = Tfinal - (Tsteps-1) * dt;
        
        U0 = setInitial(Nelm,elm_size,x,Xexc,uexc);
        U = U0;
        [ Amat,Pvmat,Pqmat,massMat,massMat_inv,mu_massMat ] = getAmat(Ord,Nelm,x,'R','L');

        for nt = 1:Tsteps
            if nt == Tsteps
                dt = dt_final;
            end
            U = RKn( Ord,x,Nelm,U,Amat,Pvmat,Pqmat,massMat,massMat_inv,mu_massMat,n_RK,dt,Time,'Dsp');
            Time = Time+dt;
        end
         Ue = get_exac( elm_size,Nelm,x,Ord_e,ir_e,x_e,Nelm_e,elm_size_e,Uexc );
         L2_ErrorStore(Ord+1,ir)= l2err_discrete( U,Ue,Ord,Nelm);
         Linf_ErrorStore(Ord+1,ir)= l8err_discrete( U,Ue,Ord,Nelm);
    end
    L2_OrderStore(Ord+1,:)=L2_OrderStore(Ord+1,:)+ErrorOrder(L2_ErrorStore(Ord+1,:));
    Linf_OrderStore(Ord+1,:)=Linf_OrderStore(Ord+1,:)+ErrorOrder(Linf_ErrorStore(Ord+1,:));
end

format short
disp(L2_OrderStore)
fprintf("\n")
disp(Linf_OrderStore)
format shortE
fprintf("\n")
disp(L2_ErrorStore)
fprintf("\n")
disp(Linf_ErrorStore)