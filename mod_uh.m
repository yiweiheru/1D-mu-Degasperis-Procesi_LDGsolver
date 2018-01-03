function [unew] = mod_uh(dl,dr,mu,Ord)
%MOD_UH Summary of this function goes here
%   Detailed explanation goes here
if Ord == 1
   unew = zeros(1,Ord+1);
   unew(1,:)=[mu-dl,mu+dr];
end
if Ord >= 2
    npt_quad = Ord+3;
    [qpt, qwt] = QuadLG( npt_quad );
    mu_un2 = zeros(3,1);
    for k = 1 : npt_quad
        mu_un2(:,1) = mu_un2(:,1) + qwt(k)*basis_1d( 2,qpt(k));
    end
    a0 = mu-dl;
    a2 = mu+dr;
    a1 = (mu*2 - ( a0*mu_un2(1,1)+a2*mu_un2(3,1) ))/mu_un2(2,1);
    xI = linspace(-1,1,Ord+1);
    unew=zeros(1,Ord+1);
    for i = 1:Ord+1
        unew(1,i)=[a0,a1,a2]*basis_1d(2,xI(i));
    end
end
end

