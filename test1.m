clear
close all

Tfinal = 20;
figure_at_time=[0,1,3,5,10,20];
ord_num = 4;
ir_num = 5;

n_RK  = 3;
period = 1;
CS =2;	% indicator of the initial data

P0 = zeros(CS,1);
Q0 = zeros(CS,1);
S0 = zeros(CS,1);
switch CS
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

UStore = zeros((ord_num+1)*(10*2^(ir_num))+1,ir_num,ord_num);
UexcStore = zeros((ord_num+1)*(10*2^(ir_num))+1,ir_num,ord_num);

fig = 0;
Ord = ord_num;
ir = ir_num;

Nelm = 10*2^(ir-1);
dx = period/Nelm;
x = 0:dx:period;
elm_size = Ord+1;

cfl = 0.1;
dt = cfl * dx;
Tsteps = floor((Tfinal-0.1*dt)/dt)+1;
dt_final = Tfinal - (Tsteps-1) * dt;

U0 = setInitial_shock(Nelm,elm_size,x,CS,period,P0,Q0,S0); % shock peakons

[ Amat,Pvmat,massMat,massMat_inv,mu_massMat ] = getAmat( Ord,Nelm,x );
U = U0;

Time = 0;
for nt = 1:Tsteps
    
    if Time == 0
        [x_to_plot,uh_to_plot]=plot_uh( U,Ord,Nelm,x);
        fig = fig+1
        figure(fig)
	      plot(x_to_plot,uh_to_plot,'-b')
        grid on
        xlabel('x')
        ylabel('u')
        Title_str = strcat('t=',num2str(Time));
        title(Title_str);
    end
    
    if nt == Tsteps-1
        dt = dt_final;
    end
    time = Time+dt;
    for p=1:size(figure_at_time,2)
    	if figure_at_time(p)-dt/2 < Time && Time <= figure_at_time(p)+dt/2
	   [x_to_plot,uh_to_plot]=plot_uh( U,Ord,Nelm,x);
	   fig = fig+1
	   figure(fig)
	   plot(x_to_plot,uh_to_plot,'or')

	   [ Uexc,~,~ ] = shock_pkns_solu( P0,Q0,S0,Nelm,elm_size, x ,period,Time,CS );
	   [x_to_plot,uh_to_plot]=plot_uh( U,Ord,Nelm,x);
	   fig = fig+1
	   figure(fig)
	   plot(x_to_plot,uh_to_plot,'-b')
	   grid on
	   xlabel('x')
	   ylabel('u')
	   Title_str = strcat('t=',num2str(Time));
	   title(Title_str);
        end
    end
end
