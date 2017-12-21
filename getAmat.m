function [ Amat,Pvmat,massMat,massMat_inv,mu_massMat ] = getAmat( Ord,Nelm,x )

% 1 represent Left, 2 represent Right

elm_size=Ord+1;

[ Mmat,Pmat,mu_Mmat,FQmat,FVmat,M_inv ] = Mat_sys( Ord,x,Nelm );

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

%compute the matrix formulate by the numerical flux
isp3=zeros(2*2*Nelm*elm_size*elm_size+1,1);
jsp3=zeros(2*2*Nelm*elm_size*elm_size+1,1);
sdat3=zeros(2*2*Nelm*elm_size*elm_size+1,1);

isp4=zeros(2*2*Nelm*elm_size*elm_size+1,1);
jsp4=zeros(2*2*Nelm*elm_size*elm_size+1,1);
sdat4=zeros(2*2*Nelm*elm_size*elm_size+1,1);

it=0;jt=0;
for ne=1:Nelm
    for i=1:elm_size
        for j=1:elm_size

            it=it+1;
            isp3(it)=(ne-1)*elm_size+i;
            jsp3(it)=(ne-1)*elm_size+j;
            sdat3(it)=0.5*FQmat(i,j,1,1,ne);
            it=it+1;
            isp3(it)=(ne-1)*elm_size+i;
            jsp3(it)=(nei(ne,1)-1)*elm_size+j;
            sdat3(it)=0.5*FQmat(i,j,1,2,ne);
            
            it=it+1;
            isp3(it)=(ne-1)*elm_size+i; 
            jsp3(it)=(nei(ne,2)-1)*elm_size+j;
            sdat3(it)=-0.5*FQmat(i,j,2,1,ne); 
            it=it+1;
            isp3(it)=(ne-1)*elm_size+i; 
            jsp3(it)=(ne-1)*elm_size+j;
            sdat3(it)=-0.5*FQmat(i,j,2,2,ne); 
            
            %rh_hat take minus
            jt=jt+1;
            isp4(jt)=(ne-1)*elm_size+i;
            jsp4(jt)=(ne-1)*elm_size+j;
            sdat4(jt)=0.5*FVmat(i,j,1,1,ne);
            jt=jt+1;
            isp4(jt)=(ne-1)*elm_size+i;
            jsp4(jt)=(nei(ne,1)-1)*elm_size+j;
            sdat4(jt)=0.5*FVmat(i,j,1,2,ne);
            
            jt=jt+1;
            isp4(jt)=(ne-1)*elm_size+i; 
            jsp4(jt)=(nei(ne,2)-1)*elm_size+j;
            sdat4(jt)=-0.5*FVmat(i,j,2,1,ne); 
            jt=jt+1;
            isp4(jt)=(ne-1)*elm_size+i; 
            jsp4(jt)=(ne-1)*elm_size+j;
            sdat4(jt)=-0.5*FVmat(i,j,2,2,ne); 
            
        end
    end
end

it=it+1;
isp3(it)=Nelm*elm_size;
jsp3(it)=Nelm*elm_size;
sdat3(it)=0;

jt=jt+1;
isp4(jt)=Nelm*elm_size;
jsp4(jt)=Nelm*elm_size;
sdat4(jt)=0;

Pqmat = PriMat + sparse(isp3,jsp3,sdat3); % r-ux=0 ,flux uh take plus
Pvmat = PriMat + sparse(isp4,jsp4,sdat4); % u-rx=q ,flux rh take minus

Amat =  mu_massMat - Pqmat * massMat_inv * Pvmat;
% Amat = - Pqmat * massMat_inv * Pvmat;


end

