function [ Uexc,P,Q ] = multi_pkns_solu( P0,Q0,Nelm,elm_size, x ,period,Time,CS )

Uexc = zeros(Nelm*elm_size,1);
P = 0;
Q = 0;


dt=1e-3;
Nt=floor((Time-0.1*dt)/dt)+1;
dt_final=Time-dt*(Nt-1);

xi=zeros(CS,CS);
A=zeros(CS,CS);
B=zeros(CS,CS);

% ODE solver with Euler method
for it=1:Nt
    
    if it == Nt-1
        dt = dt_final;
    end
    
    for i=1:CS
        for j=1:CS
            xi( i,j ) = Q0(i)-Q0(j);
            xi( i,j ) = xi(i,j)-floor(xi(i,j)/period)*period;
            % Generate the matrix A for g
            A(i,j) = 0.5*(xi(i,j)-1/2 )^2 + 23/24;
            % Generate the matrix B for g'
            if xi(i,j) == 0 || xi(i,j) == period
                B(i,j) = 0;
            else
                B(i,j) = -2*P0(i) * (xi(i,j)-1/2); % The G-DP equation should change the parameter (-1 --> -2)
            end
        end
    end
    
    Q0 = Q0+A*P0.*dt;
    P0 = P0+B*P0.*dt;
end

Q=Q0;
P=P0;

for ne = 1:Nelm
    for i = 1:elm_size
        xtemp = x(ne)+(x(ne+1)-x(ne))*(i-1)/(elm_size-1);
        num = (ne-1)*elm_size+i;
        switch CS
            case 1
                x_tilde1 = xtemp - Q(1);
                x_hat1 = x_tilde1 - floor(x_tilde1/period)*period ;
                psi1 = P(1) * (0.5 * ( x_hat1-1/2 )^2 + 23/24);
                
                Uexc(num) = psi1 ;
            case 2
                x_tilde1 = xtemp - Q(1);
                x_hat1 = x_tilde1 - floor(x_tilde1/period)*period ;
                psi1 = P(1) * (0.5 * ( x_hat1-1/2 )^2 + 23/24);
                
                x_tilde2 = xtemp - Q(2);
                x_hat2 = x_tilde2 - floor(x_tilde2/period)*period ;
                psi2 = P(2) *( 0.5 * (x_hat2-1/2 )^2 + 23/24);
                
                Uexc(num) = psi1 + psi2 ;
            case 3
                x_tilde1 = xtemp - Q(1);
                x_hat1 = x_tilde1 - floor(x_tilde1/period)*period ;
                psi1 = P(1) * (0.5 * ( x_hat1-1/2 )^2 + 23/24);
                
                x_tilde2 = xtemp - Q(2);
                x_hat2 = x_tilde2 - floor(x_tilde2/period)*period ;
                psi2 = P(2) *( 0.5 * (x_hat2-1/2 )^2 + 23/24);
                
                x_tilde3 = xtemp - Q(3);
                x_hat3 = x_tilde3 - floor(x_tilde3/period)*period ;
                psi3 = P(3) *( 0.5 * (x_hat3-1/2 )^2 + 23/24);
                
                Uexc(num) = psi1 + psi2 + psi3;
        end
    end
end
end

