function [ U0 ]  =  setInitial_shock(Nelm,elm_size,x,CS,period,P0,Q0,S0)

% CS is the indicator of the type of initial data U0



U0      = zeros(Nelm*elm_size,1);

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
            
            x_tilde(k) = xtemp - Q0(k);
            
            x_hat(k)   = x_tilde(k) - floor(x_tilde(k)/period)*period;
            
            % descriptions:
            
            % g(x)  = 1/2*(x-1/2)^2+23/24;
            
            % g'(x) = x-1/2, 0<x<1 | g''(x) = 1, 0<x<1
            
            %         0    , x=0   |          0, x=0
            
            if x_hat(k) == 0 || x_hat(k) == period
                
                psi(k) = P0(k)*(0.5*(x_hat(k)-1/2)^2 + 23/24);
                
            else
                
                psi(k) = P0(k)*(0.5*(x_hat(k)-1/2)^2 + 23/24) + S0(k)*(x_hat(k)-1/2);
                
            end
            
        end
        
        U0(num) = sum(psi);
        
    end
    
end

end
