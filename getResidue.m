function [ Residue ] = getResidue( Ord,x,Nelm,U,Amat,Pvmat,massMat,massMat_inv,mu_massMat,Time )

elm_size=Ord+1;

V = Amat\( (massMat)*U );

Q = -1*massMat_inv*Pvmat*V ;

uh = uhTransform( Nelm,elm_size,U );
qh = uhTransform( Nelm,elm_size,Q );

Residue = residue( x,Nelm,Ord,uh,qh,Time );

end

