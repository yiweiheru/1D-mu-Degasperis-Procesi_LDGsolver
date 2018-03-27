function [ Amat,Pvmat,Pqmat,massMat,massMat_inv,mu_massMat ] = getAmat( Ord,Nelm,x,flux_q,flux_v )

% 1 represent Left, 2 represent Right

elm_size=Ord+1;

[ Mmat,Pmat,mu_Mmat,Fmat,M_inv ] = Mat_sys( Ord,x,Nelm );

% compute the inverse of mass matrix
isp=zeros(Nelm*elm_size*elm_size+1,1);
jsp=zeros(Nelm*elm_size*elm_size+1,1);
sdat=zeros(Nelm*elm_size*elm_size+1,1);

it=0;
for ne=1:Nelm
    Jacob=(x(ne+1)-x(ne))/2;
    for i=1:elm_size
        for j=1:elm_size
            it=it+1;
            isp(it)=(ne-1)*elm_size+i;
            jsp(it)=(ne-1)*elm_size+j;
            sdat(it)=M_inv(i,j)/Jacob;
        end
    end
end
it=it+1;
isp(it)=Nelm*elm_size;
jsp(it)=Nelm*elm_size;
sdat(it)=0;

massMat_inv=sparse(isp,jsp,sdat);

%compute the mu_mass matrix
isp=zeros(Nelm*Nelm*elm_size*elm_size+1,1);
jsp=zeros(Nelm*Nelm*elm_size*elm_size+1,1);
sdat=zeros(Nelm*Nelm*elm_size*elm_size+1,1);

it=0;
for neI=1:Nelm
    for neJ=1:Nelm
        for i=1:elm_size
            for j=1:elm_size
                it=it+1;
                isp(it)=(neI-1)*elm_size+i;
                jsp(it)=(neJ-1)*elm_size+j;
                sdat(it)=mu_Mmat(i,j,neI,neJ);
            end
        end
    end
end
it=it+1;
isp(it)=Nelm*elm_size;
jsp(it)=Nelm*elm_size;
sdat(it)=0;

mu_massMat=sparse(isp,jsp,sdat);

%compute the mass matrix & prime matrix

%compute the mass matrix & prime matrix
isp1=zeros(Nelm*elm_size*elm_size+1,1);
jsp1=zeros(Nelm*elm_size*elm_size+1,1);
sdat1=zeros(Nelm*elm_size*elm_size+1,1);

isp2=zeros(Nelm*elm_size*elm_size+1,1);
jsp2=zeros(Nelm*elm_size*elm_size+1,1);
sdat2=zeros(Nelm*elm_size*elm_size+1,1);


it=0;jt=0;
for ne=1:Nelm
    for i=1:elm_size
        for j=1:elm_size
            it=it+1;
            isp1(it)=(ne-1)*elm_size+i;
            jsp1(it)=(ne-1)*elm_size+j;
            sdat1(it)=Mmat(i,j,ne);
            
            jt=jt+1;
            isp2(jt)=(ne-1)*elm_size+i;
            jsp2(jt)=(ne-1)*elm_size+j;
            sdat2(jt)=Pmat(i,j,ne);
        end
    end
end

it=it+1;
isp1(it)=Nelm*elm_size;
jsp1(it)=Nelm*elm_size;
sdat1(it)=0;

jt=jt+1;
isp2(jt)=Nelm*elm_size;
jsp2(jt)=Nelm*elm_size;
sdat2(jt)=0;

massMat=sparse(isp1,jsp1,sdat1);
PriMat=sparse(isp2,jsp2,sdat2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Pqmat = PriMat + getAmat_flux(Nelm,elm_size,Fmat,flux_q);
Pvmat = PriMat + getAmat_flux(Nelm,elm_size,Fmat,flux_v);

Amat =  mu_massMat - Pqmat * massMat_inv * Pvmat;

% % Amat*Qh = Uh;
% Amat = massMat_inv * (-mu_massMat*massMat_inv*Pvmat + Pqmat);

end

