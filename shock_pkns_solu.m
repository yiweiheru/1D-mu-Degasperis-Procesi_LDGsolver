function [ Uexc,P,Q,S ] = shock_pkns_solu( P0,Q0,S0,Nelm,elm_size, x ,period,Time,CS )

Lam = 3;
    
    dt=1e-3;
    Nt=floor((Time-0.1*dt)/dt)+1;
    dt_final=Time-dt*(Nt-1);
    
    xi = zeros(CS,CS);
    G0 = zeros(CS,CS);
    G1 = zeros(CS,CS);
	G2 = ones(CS,CS); 
    
    % ODE solver with Euler method
    for it=1:Nt
        
        if it == Nt
            dt = dt_final;
        end
        
        for i=1:CS
            for j=1:CS
                xi( i,j ) = Q0(i)-Q0(j);
                xi( i,j ) = xi(i,j)-floor(xi(i,j)/period)*period;
                % Generate the matrix A for g
                G0(i,j) = 0.5*(xi(i,j)-1/2 )^2 + 23/24;
                % Generate the matrix B for g'
                if xi(i,j) == 0 || xi(i,j) == period
                    G1(i,j) = 0;
                else
                    G1(i,j) = (xi(i,j)-1/2); 
                end
            end
        end
        
        Q0 = Q0 + (G0*P0 + G1*S0)*dt;
        P0 = P0 + (Lam-1)*(S0.*(G2*P0) - P0.*(G1*P0 + G2*S0))*dt;
		S0 = S0 + -(Lam-2)*(S0.*(G1*P0 + G2*S0))*dt;
    end
    Q = Q0;
    P = P0;
	S = S0;
	
	
    Uexc    = zeros(Nelm*elm_size,1);
	x_tilde = zeros(1,CS);
	x_hat   = zeros(1,CS);
	psi     = zeros(1,CS);
    for ne = 1:Nelm
        for i = 1:elm_size
            if elm_size == 1
                xtemp = x(ne);
            else
                xtemp = x(ne)+(x(ne+1)-x(ne))*(i-1)/(elm_size-1);
            end
            num = (ne-1)*elm_size+i;
            for k = 1:CS
				x_tilde(k) = xtemp - Q(k);
				x_hat(k)   = x_tilde(k) - floor(x_tilde(k)/period)*period;
				if (x_hat(k) == 0 || x_hat(k) == period)
					psi(k) = P(k)*(0.5*(x_hat(k)-1/2)^2 + 23/24);
				else
					psi(k) = P(k)*(0.5*(x_hat(k)-1/2)^2 + 23/24) + S(k)*(x_hat(k)-1/2);
                end
            end
            Uexc(num) = sum(psi);
        end
    end
end
