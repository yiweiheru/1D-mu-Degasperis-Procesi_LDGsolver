function [ U0 ]  =  setInitial_shock(Nelm,elm_size,x,CS,period,P0,Q0,S0)
% CS is the indicator of the type of initial data U0

U0    = zeros(Nelm*elm_size,1);
x_tilde = zeros(1,CS);
x_hat   = zeros(1,CS);
psi     = zeros(1,CS);

for ne = 1:Nelm
    for i = 1:elm_size
        xtemp = x(ne)+(x(ne+1)-x(ne))*(i-1)/(elm_size-1);
        num = (ne-1)*elm_size+i;
        
        for k = 1:CS
            x_tilde(k) = xtemp - Q0(k);
            x_hat(k)   = x_tilde(k) - floor(x_tilde(k)/period)*period;
            if x_hat(k) == 0 || x_hat(k) == period
                psi(k) = P0(k)*(0.5*(x_hat(k)-1/2)^2 + 23/24);
            else
                psi(k) = P0(k)*(0.5*(x_hat(k)-1/2)^2 + 23/24) + S0(k)*(x_hat(k)-1/2);
            end
            U0(num) = sum(psi);
        end
    end
end
end




