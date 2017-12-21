function expr = basis_1d( ord,var )

expr=zeros(ord+1,1);

switch (ord)
    case 1 % Linear polynomial 
        expr(1) = (1-var)/2;
        expr(2) = (1+var)/2;   
    case 2 % Quadratic polynomial    
        expr(1) = -1/2*var*(1-var);
        expr(2) = (1-var)*(1+var);
        expr(3) = 1/2*var*(1+var);        
    case 3 % Third order polynomial
        expr(1) = -9/16*(1-var)*(1/3+var)*(1/3-var);
        expr(2) = 27/16*(1-var)*(1+var)*(1/3-var);
        expr(3) = 27/16*(1-var)*(1+var)*(1/3+var);
        expr(4) = -9/16*(1+var)*(1/3+var)*(1/3-var);
    case 4 % Fourth order polynomial
        expr(1) = (2*var-1)*(var-1)*var*(2*var+1)/6;
        expr(2) = -4*(2*var-1)*(var-1)*var*(var+1)/3;
        expr(3) = (2*var-1)*(var-1)*(var+1)*(2*var+1);
        expr(4) = -4*(var-1)*var*(var+1)*(2*var+1)/3;
        expr(5) = (2*var-1)*var*(var+1)*(2*var+1)/6;
    case 5 % Fifth order polynomial
        expr(1) = -(5*var+3)*(5*var+1)*(5*var-1)*(5*var-3)*(var-1)/768;
        expr(2) = 25*(5*var+1)*(var+1)*(5*var-1)*(5*var-3)*(var-1)/768;
        expr(3) = -25*(5*var+3)*(var+1)*(5*var-1)*(5*var-3)*(var-1)/384;
        expr(4) = 25*(5*var+3)*(5*var+1)*(var+1)*(5*var-3)*(var-1)/384;
        expr(5) = -25*(5*var+3)*(5*var+1)*(var+1)*(5*var-1)*(var-1)/768;
        expr(6) = (5*var+3)*(5*var+1)*(var+1)*(5*var-1)*(5*var-3)/768;
end % Switch on order
% switch (ord)
%     case 1 % Linear polynomial 
%         expr(1) = 1;
%         expr(2) = var/2;   
%     case 2 % Quadratic polynomial    
%         expr(1) =1;
%         expr(2) = var/2;
%         expr(3) = (var/2)^2;        
% end % Switch on order

return
