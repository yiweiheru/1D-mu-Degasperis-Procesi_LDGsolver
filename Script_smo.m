
% run simu_exac.m
close all
clear

load exact.mat
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
        
        Time = 0;
        Nelm = 10*2^(ir-1)+1;
        dx = period/Nelm;
        x = 0:dx:period;
        elm_size = Ord+1;
        
        dt = cfl * dx;
        Tsteps = floor((Tfinal-0.1*dt)/dt)+1;
        dt_final = Tfinal - (Tsteps-1) * dt;
        
        U0 = setInitial_smo(Nelm,elm_size,x);
        U = U0;
        [ Amat,Pvmat,Pqmat,massMat,massMat_inv,mu_massMat ] = getAmat(Ord,Nelm,x,'C','C');
      
        for nt = 1:Tsteps
            if nt == Tsteps
                dt = dt_final;
            end
            U = RKn( Ord,x,Nelm,U,Amat,Pvmat,Pqmat,massMat,massMat_inv,mu_massMat,n_RK,dt,Time,'Csv');
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
disp(L2_ErrorStore)
fprintf("\n")
disp(Linf_ErrorStore)

fprintf("rank of Pq*massMat_inv*Pv is %i \n", rank(full(mu_massMat - Amat)));
fprintf("rank of Amat:%i, size of Amat:%i \n", rank(full(Amat)),size(full(Amat),1));
fprintf("condition number of Amat is %e \n",  condest(full(Amat)));

