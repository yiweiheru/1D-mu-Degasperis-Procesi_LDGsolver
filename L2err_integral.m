function [ L2err ] = L2err_integral( Uh,Ue,Ord,Nelm,x )

elm_size=Ord+1;
U_diff=Ue-Uh;
u=uhTransform( Nelm,elm_size,U_diff );

npt_quad=Ord+2;
[qpt, qwt] = QuadLG(npt_quad);

un=zeros(elm_size,npt_quad);
for k = 1 : npt_quad
    un(:,k)=basis_1d(Ord,qpt(k));
end

Val=0;
for ne=1:Nelm
    Jaco=(x(ne+1)-x(ne))/2;
    for ik=1:npt_quad
        Val=Val+(u(ne,:)*un(:,ik))^2*Jaco*qwt(ik);
    end
end

L2err=sqrt(Val);

end
