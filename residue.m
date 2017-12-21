function [ Residue ] = residue( x,Nelm,Ord,uh,qh,~ )

elm_size = Ord+1;

npt_quad = Ord+3;
[qpt,qwt] = QuadLG(npt_quad);

un     = zeros( elm_size,npt_quad );
un_der = zeros( elm_size,npt_quad );
uf     = zeros( elm_size,2 );
mu_un  = zeros(elm_size,1);

for k = 1 : npt_quad
    un(:,k)     = basis_1d( Ord,qpt(k) );
    un_der(:,k) = basisDer_1d( Ord,qpt(k) );
    mu_un(:,1)  = mu_un(:,1) + qwt(k)*un(:,k);
end

uf(:,1) = basis_1d(Ord,-1); %left boundary of elment
uf(:,2) = basis_1d(Ord,1);  %right boundary of element

Residue   = zeros(Nelm*elm_size,1);
Intergral = zeros(Nelm*elm_size,1);
Flux_L    = zeros(Nelm*elm_size,1);
Flux_R    = zeros(Nelm*elm_size,1);

%when implement the average of un, Jacobbi determinent is needed
mu_uh = 0;
for ne = 1 : Nelm
    Jaco = (x(ne+1)-x(ne))/2;
    for j = 1 : elm_size
        mu_uh = mu_uh+mu_un(j,1)*uh(ne,j)*Jaco;
    end
end


%imply the periodic boundary condition to get the neighbourhood.
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




for ne=1:Nelm
    Jaco = (x(ne+1)-x(ne))/2;
    for ik=1:npt_quad
        u=0; q=0;
        for j=1:elm_size
            u=u+un(j,ik)*uh(ne,j);
            q=q+un(j,ik)*qh(ne,j);
        end
        f=0.5*u^2;
        for i=1:elm_size
            num=(ne-1)*elm_size+i;
            Intergral(num)=Intergral(num)+qwt(ik)*(un_der(i,ik)*(-1*f) + un(i,ik)*3*mu_uh*q*Jaco);
        end
    end
end

% get the maxmimum of u
uh_max = max(max(abs(uh)));

for ne=1:Nelm
    
    uLp = uh(ne,:)*uf(:,1);
    uLm = uh(nei(ne,1),:)*uf(:,2);
    uRp = uh(nei(ne,2),:)*uf(:,1);
    uRm = uh(ne,:)*uf(:,2);
    
    
    fhat_L = 0.5*(0.5*uLp^2 + 0.5*uLm^2 - uh_max*(uLp - uLm));
    fhat_R = 0.5*(0.5*uRp^2 + 0.5*uRm^2 - uh_max*(uRp - uRm));
    
    
    for i=1:elm_size
        num=(ne-1)*elm_size+i;
        Flux_L(num)=Flux_L(num)+uf(i,1)*(-1*fhat_L);
        Flux_R(num)=Flux_R(num)+uf(i,2)*(-1*fhat_R);
    end
end

Residue=Residue-(Intergral-Flux_R+Flux_L);

end


