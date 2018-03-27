function [ U0 ]  =  setInitial_peak(Nelm,elm_size,x,CS,period,P0,Q0)
% CS is the indicator of the type of initial data U0



U0 = zeros(Nelm*elm_size,1);
for ne = 1:Nelm
    for i = 1:elm_size
        xtemp = x(ne)+(x(ne+1)-x(ne))*(i-1)/(elm_size-1);
        num = (ne-1)*elm_size+i;
        
        switch CS
            %-------------------------One peakon---------------------------
            case 1
                q1 = Q0(1);p1 = P0(1);
                
                x_tilde = xtemp - q1 ;
                x_hat   = x_tilde - floor(x_tilde/period) * period ;
                psi     = p1 * (0.5 *  ( x_hat-1/2 )^2 + 23/24);
                U0(num) = psi;
%                 x_hat = (x_tilde-(-period/2))-floor((x_tilde-(-period/2))/period)*period+(-period/2);
                
%                 U0(num) = sin(x_hat);
%                 U0(num) = p1/26*(12*x_hat^2+23);
%                 % scaling!
%                 argu_nv =p1*period/(1+1/12*period^3);
%                 argu_m = p1-period^2/8*argu_nv;
%                 x_tilde = xtemp - q1 + 0.5*period;
%                 x_hat = (x_tilde-(-0.5*period))-floor((x_tilde-(-0.5*period))/period)*period+(-0.5*period); % Periodically in R
%                 U0(num) = argu_m+(argu_nv/2)*x_hat^2;
                %----------------------Two peakons-------------------------
            case 2
                q1 = Q0(1);p1 = P0(1);
                q2 = Q0(2);p2 = P0(2);
                
                x_tilde1 = xtemp - q1;
                x_hat1 = x_tilde1 - floor(x_tilde1/period) * period ;
                psi1 = p1 * (0.5 *  ( x_hat1-1/2 )^2 + 23/24);
                
                x_tilde2 = xtemp - q2;
                x_hat2 = x_tilde2 - floor(x_tilde2/period) * period ;
                psi2 = p2 *( 0.5 *  (x_hat2-1/2 )^2 + 23/24);
                
                U0(num) = psi1 + psi2;
                %----------------------Three peakons-----------------------
            case 3
                q1 = Q0(1);p1 = P0(1);
                q2 = Q0(2);p2 = P0(2);
                q3 = Q0(3);p3 = P0(3);
                
                x_tilde1 = xtemp - q1;
                x_hat1 = x_tilde1 - floor(x_tilde1/period)*period ;
                psi1 = p1 * (0.5 *  ( x_hat1-1/2 )^2 + 23/24);
                
                x_tilde2 = xtemp - q2;
                x_hat2 = x_tilde2 - floor(x_tilde2/period)*period ;
                psi2 = p2 *( 0.5 *  (x_hat2-1/2 )^2 + 23/24);
                
                x_tilde3 = xtemp - q3;
                x_hat3 = x_tilde3 - floor(x_tilde3/period)*period ;
                psi3 = p3 *( 0.5 *  (x_hat3-1/2 )^2 + 23/24);
                
                U0(num) = psi1 + psi2 +psi3;
                
                
        end
    end
end
end


