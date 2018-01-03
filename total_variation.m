function [TV] = total_variation(U,Nelm,Ord,x)
%TOTAL_VARIATION Summary of this function goes here
%   Detailed explanation goes here
elm_size = Ord+1;
uh = uhTransform( Nelm,elm_size,U );

npt_quad  = Ord+3;
[qpt,qwt] = QuadLG(npt_quad);

mu_un  = zeros( elm_size,1);

for k = 1 : npt_quad
    mu_un(:,1) = mu_un(:,1) + qwt(k)*basis_1d( Ord,qpt(k) );
end

mu_uh = zeros(Nelm);
for ne = 1:Nelm
    Jaco = (x(ne+1)-x(ne))/2;
    mu_uh(ne) = (1/(2*Jaco))*uh(ne,:)*mu_un*Jaco;
end

nei=zeros(Nelm,2);
for ne=1:Nelm
    if ne==1
        nei(ne,1)=Nelm;
        nei(ne,2)=2;
    elseif ne==Nelm
        nei(ne,1)=Nelm-1;
        nei(ne,2)=1;
    else
        nei(ne,1)=ne-1;
        nei(ne,2)=ne+1;
    end
end

TV = 0;
for ne = 1:Nelm
    TV = TV + abs(mu_uh(nei(ne,2)) - mu_uh(ne));
end

end

