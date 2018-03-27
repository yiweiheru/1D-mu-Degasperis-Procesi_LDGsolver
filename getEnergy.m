function [ H ] = getEnergy( Uh,Ord,Nelm,x)

elm_size=Ord+1;
U=Uh;
u=uhTransform( Nelm,elm_size,U );

npt_quad=Ord+2;
[qpt, qwt] = QuadLG(npt_quad);
un=zeros(elm_size,npt_quad);
mu_un=zeros(elm_size,1);
for ik = 1 : npt_quad
    un(:,ik)=basis_1d(Ord,qpt(ik));
    mu_un(:,1)=mu_un(:,1)+qwt(ik)*un(:,ik);
end
%--------------------------compute energy of r-----------------------------
Val=0;
for ne=1:Nelm
    Jaco=(x(ne+1)-x(ne))/2;
    for ik=1:npt_quad
        Val=Val+(u(ne,:)*un(:,ik))^2*Jaco*qwt(ik);
    end
end
% -------------------------compute mean of u-------------------------------
mean_u=0;
for ne=1:Nelm 
    Jaco=(x(ne+1)-x(ne))/2;
    mean_u = mean_u + u(ne,:)*mu_un(:,1)*Jaco;
end
%--------------------------------------------------------------------------
H0 = mean_u;
H1 = (Val);
H = [H0,H1];
end

