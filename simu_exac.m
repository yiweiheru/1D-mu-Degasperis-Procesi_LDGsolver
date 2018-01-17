% In this script, we utilize the sufficient small mesh grid to approximate
% the exact solution for specific initial condition.
clear
Ord_e = 3;
ir_e     = 6;
cfl     = 0.1;
n_RK    = 3;
period  = 1;
Tfinal=0.5;

Time = 0;
Nelm = 10*2^(ir_e-1);
dx = period/Nelm;
x = 0:dx:period;
elm_size = Ord_e+1;

dt = cfl * dx;
Tsteps  = floor((Tfinal-0.1*dt)/dt)+1;
dt_final = Tfinal - (Tsteps-1) * dt;

U0 = setInitial_smo(Nelm,elm_size,x);
U = U0;
[ Amat,Pvmat,massMat,massMat_inv,mu_massMat ] = getAmat(Ord_e,Nelm,x);

for nt = 1:Tsteps
    if nt == Tsteps-1
        dt = dt_final;
    end
    
    U = RKn( Ord_e,x,Nelm,U,Amat,Pvmat,massMat,massMat_inv,mu_massMat,n_RK,dt,Time);
    Time = Time+dt;
end

Uexc = U;

save('exact.mat','Ord_e','ir_e','Uexc','period','Tfinal','cfl')
% plot_uh( Uexc,Ord_e,(10*2^(ir_e-1)),x ,"numerical")