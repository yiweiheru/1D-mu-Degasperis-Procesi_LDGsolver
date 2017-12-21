function [ Mmat,Pmat,mu_Mmat,FQmat,FVmat,M_inv ] = Mat_sys( Ord,x,Nelm )
%    In this projection of DP equation, we choose a simple way
%    to handle the mass matrix via using uniform mesh.

%    1 represent Left, 2 represent Right

elm_size = Ord+1;
npt_quad = Ord+2;
[qpt, qwt] = QuadLG( npt_quad );

un = zeros( elm_size,npt_quad );
un_der = zeros( elm_size,npt_quad );
mu_un = zeros( elm_size,1 );
for k=1 : npt_quad
    un(:,k)        = basis_1d( Ord,qpt(k) );
    un_der(:,k) = basisDer_1d( Ord,qpt(k) );
    mu_un(:,1) = mu_un(:,1) + qwt(k)*un(:,k);
end

Mmat = zeros(elm_size,elm_size,Nelm);
Pmat = zeros(elm_size,elm_size,Nelm);
mu_Mmat = zeros(elm_size,elm_size,Nelm,Nelm);

for ne=1:Nelm
    Joca=(x(ne+1)-x(ne))/2;
    for k=1:npt_quad
        for j=1:elm_size
            for i=1:elm_size
                Mmat(i,j,ne)=Mmat(i,j,ne)+qwt(k)*un(i,k)*un(j,k)*Joca;
                Pmat(i,j,ne)=Pmat(i,j,ne)+qwt(k)*un_der(i,k)*un(j,k);
            end
        end
    end
end

%to avoid the bad condition number,we compute the inv_massMat
M0mat=zeros(elm_size,elm_size);
for k=1:npt_quad
    for j=1:elm_size
        for i=1:elm_size
            M0mat(i,j)=M0mat(i,j)+qwt(k)*un(i,k)*un(j,k);
        end
    end
end
M_inv=inv(M0mat);

for neI=1:Nelm
    Joca1=(x(neI+1)-x(neI))/2;
    for neJ=1:Nelm
        Joca2=(x(neJ+1)-x(neJ))/2;
        for j=1:elm_size
            for i=1:elm_size
                mu_Mmat(i,j,neI,neJ)=mu_Mmat(i,j,neI,neJ)+mu_un(i,1)*Joca1*mu_un(j,1)*Joca2;
            end
        end
    end
end
% %imply the periodic boundary condition to get the neighbourhood.
% nei=zeros(Nelm,2);
% for ne=1:Nelm
%     if ne==1
%         nei(ne,1)=  Nelm;
%         nei(ne,2)=2;
%     elseif ne==Nelm
%         nei(ne,1)=Nelm-1;
%         nei(ne,2)=1;
%     else
%         nei(ne,1)=ne-1;
%         nei(ne,2)=ne+1;
%     end
% end

% Matrix fomulated by the numerical flux, they are different from the
% choice of the flux

uf=zeros(elm_size,2);
uf(:,1)=basis_1d(Ord,-1);
uf(:,2)=basis_1d(Ord,1);

FQmat = zeros( elm_size,elm_size,2,2,Nelm );   % uh takes plus
FVmat = zeros( elm_size,elm_size,2,2,Nelm );   % rh takes minus

for ne = 1:Nelm
    for j = 1:elm_size
        for i = 1:elm_size
            FQmat(i,j,1,1,ne) = uf(i,1)*uf(j,1);
            FQmat(i,j,1,2,ne) = uf(i,1)*uf(j,2);
            FQmat(i,j,2,1,ne) = uf(i,2)*uf(j,1);
            FQmat(i,j,2,2,ne) = uf(i,2)*uf(j,2);
            
            FVmat(i,j,1,1,ne) = uf(i,1)*uf(j,1);
            FVmat(i,j,1,2,ne) = uf(i,1)*uf(j,2);
            FVmat(i,j,2,1,ne) = uf(i,2)*uf(j,1);
            FVmat(i,j,2,2,ne) = uf(i,2)*uf(j,2);
        end
    end
end

end

