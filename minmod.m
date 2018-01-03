function [v] = minmod(v1,v2,v3,h,M)
%MINMOD Summary of this function goes here
%   Detailed explanation goes here

if abs(v1) <= M*h^2
    v = v1;
elseif (sign(v1) == sign(v2)) && (sign(v1) == sign(v3))
    v = sign(v1)*min(abs([v1,v2,v3]));
else
    v = 0;
end
end

