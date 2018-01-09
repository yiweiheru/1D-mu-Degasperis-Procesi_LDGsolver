function [Umod] = TVB_limiter( U,Ord,Nelm,x,P0 )
% use the TVB limiter to remove the oscilation near the 
% discontinuouity of u_h

% M=2/3 * max(abs(P0));
M = 50;
elm_size = Ord + 1;
uh = uhTransform( Nelm,elm_size,U );

npt_quad  = Ord+3;
[qpt,qwt] = QuadLG(npt_quad);
un     = zeros( elm_size,2 );
mu_un  = zeros( elm_size,1);

un(:,1) = basis_1d(Ord,-1);
un(:,2) = basis_1d(Ord, 1);
for k = 1 : npt_quad
    mu_un(:,1) = mu_un(:,1) + qwt(k)*basis_1d( Ord,qpt(k) );
end

mu_uh = zeros(Nelm);
df_L  = zeros(Nelm);
df_R  = zeros(Nelm);
for ne = 1:Nelm
    Jaco = (x(ne+1)-x(ne))/2;
    mu_uh(ne) = (1/(2*Jaco))*uh(ne,:)*mu_un*Jaco;
    df_R(ne)  = uh(ne,:)*un(:,2) - mu_uh(ne);
    df_L(ne)  = mu_uh(ne) - uh(ne,:)*un(:,1);
end

% the core of the TVB limiter:minmod function
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

umod = zeros(Nelm,elm_size);
for ne = 1:Nelm
    h = x(ne+1)-x(ne);
    dl = minmod(df_L(ne),mu_uh(nei(ne,2))-mu_uh(ne),mu_uh(ne)-mu_uh(nei(ne,1)),h,M);
    dr = minmod(df_R(ne),mu_uh(nei(ne,2))-mu_uh(ne),mu_uh(ne)-mu_uh(nei(ne,1)),h,M);
    
    if (dl == df_L(ne)) && (dr == df_R(ne))
        umod(ne,:) = uh(ne,:);
    else
        umod(ne,:) = mod_uh(dl,dr,mu_uh(ne),Ord);
    end
end
Umod=zeros(Nelm*elm_size,1);
for ne = 1:Nelm
    for i = 1:elm_size
        num = (ne-1)*elm_size+i;
        Umod(num)= umod(ne,i);
    end
end



end