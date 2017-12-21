function expr = basisDer_1d( ord,var )

expr=zeros(ord+1,1);

switch (ord)
    case 1 % For linear polynomial
        expr(1) = -1/2;
        expr(2) = 1/2;
    case 2 % For quadratic polynomial
        expr(1) = var-1/2;
        expr(2) = -2*var;
        expr(3) = var+1/2;
    case 3 % For third order polynomial
        expr(1) = (-27/16*var^2 + 9/8*var + 1/16);
        expr(2) = (81/16*var^2 - 9/8*var - 27/16);
        expr(3) = (-81/16*var^2 - 9/8*var + 27/16);
        expr(4) = (27/16*var^2 + 9/8*var - 1/16);
    case 4 % For fourth order polynomial
        expr(1) = (4*var-1)*(4*var^2 - 2*var -1)/6;
        expr(2) = -32/3*var^3 + 4*var^2 + 16/3*var -4/3;
        expr(3) = 2*var*(8*var^2-5);
        expr(4) = -32/3*var^3 - 4*var^2 + 16/3*var +4/3;
        expr(5) = (4*var+1)*(4*var^2 + 2*var -1)/6;
    case 5 % For fifth order polynomial
        expr(1) = -3125/768*var^4 + 625/192*var^3 + 125/128*var^2 - ...
            125/192*var - 3/256;
        expr(2) = 15625/768*var^4 - 625/64*var^3 - 1625/128*var^2 + ...
            325/64*var + 125/768;
        expr(3) = -15625/384*var^4 + 625/96*var^3 + 2125/64*var^2 - ...
            425/96*var - 375/128;
        expr(4) = 15625/384*var^4 + 625/96*var^3 - 2125/64*var^2 - ...
            425/96*var + 375/128;
        expr(5) = -15625/768*var^4 - 625/64*var^3 + 1625/128*var^2 + ...
            325/64*var - 125/768;
        expr(6) = 3125/768*var^4 + 625/192*var^3 - 125/128*var^2 - ...
            125/192*var + 3/256;
end
% switch (ord)
%     case 0
%         expr(1) = 0;
%     case 1 % Linear polynomial
%         expr(1) = 0;
%         expr(2) = 1/2;
%     case 2 % Quadratic polynomial
%         expr(1) = 0;
%         expr(2) = 1/2;
%         expr(3) = (var/2);
%     case 3 % Third order polynomial
%         expr(1) = 0;
%         expr(2) = 1/2;
%         expr(3) = (var/2);
%         expr(4) = 3/2*(var/2)^2;
%     case 4 % Third order polynomial
%         expr(1) = 0;
%         expr(2) = 1/2;
%         expr(3) = (var/2);
%         expr(4) = 3/2*(var/2)^2;
%         expr(5) = 2*(var/2)^3;
% end

return


