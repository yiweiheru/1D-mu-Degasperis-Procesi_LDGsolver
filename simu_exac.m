% In this script, we utilize the sufficient small mesh grid to approximate
% the exact solution for specific initial condition:
% 0.03*sin(2*pi*x)+0.5 in periodic [0,1];

clear
close all

Ord_e  = 4;
ir_e   = 5;
cfl    = 0.05;
n_RK   = 4;
period = 1;
Tfinal = 0.2;

Time   = 0;
Nelm_e = 10*2^(ir_e-1)+1;
dx     = period/Nelm_e;
x_e    = 0:dx:period;
elm_size_e = Ord_e+1;

dt = cfl * dx;
Tsteps  = floor((Tfinal-0.1*dt)/dt)+1;
dt_final = Tfinal - (Tsteps-1) * dt;

U0 = setInitial_smo(Nelm_e,elm_size_e,x_e);
U = U0;
plot_uh( U0,Ord_e,Nelm_e,x_e ,"exact")
[ Amat,Pvmat,Pqmat,massMat,massMat_inv,mu_massMat ] = getAmat(Ord_e,Nelm_e,x_e);

for nt = 1:Tsteps
    if nt == Tsteps
        dt = dt_final;
    end
    
    U = RKn( Ord_e,x_e,Nelm_e,U,Amat,Pvmat,Pqmat,massMat,massMat_inv,mu_massMat,n_RK,dt,Time);
    Time = Time+dt;
end

Uexc = U;

save('exact.mat','Ord_e','ir_e','x_e','Nelm_e','elm_size_e','Uexc','period','Tfinal','cfl')
hold on
plot_uh( Uexc,Ord_e,Nelm_e,x_e ,"numerical")