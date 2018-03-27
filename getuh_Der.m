function [ R ] = getuh_Der( Ord,x,Nelm,U,massMat_inv )

elm_size=Ord+1;
uh=uhTransform(Nelm,elm_size,U);

[ Residue1 ] = residue1( x,Nelm,Ord,uh );
R=massMat_inv*Residue1;

end

