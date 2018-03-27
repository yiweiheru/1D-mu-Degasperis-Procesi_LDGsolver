function [ Val ] = evalue_uh( Ord_e,Nelm_e,x_e,uh_exc,x0 )

nei=zeros(Nelm_e,2);
for ne=1:Nelm_e
    if ne==1
        nei(ne,1)=Nelm_e;
        nei(ne,2)=2;
    elseif ne==Nelm_e
        nei(ne,1)=Nelm_e-1;
        nei(ne,2)=1;
    else
        nei(ne,1)=ne-1;
        nei(ne,2)=ne+1;
    end
end

for ne=1:Nelm_e
    if x0==x_e(ne) 
        unr=basis_1d(Ord_e,-1);
        unl=basis_1d(Ord_e,1);
        Val=(uh_exc(ne,:)*unr+uh_exc(nei(ne,1),:)*unl)/2;
    elseif x0>x_e(ne) && x0<x_e(ne+1)
        xi=(x0-x_e(ne))/(x_e(ne+1)-x_e(ne))*2+(-1);
        un=basis_1d(Ord_e,xi);
        Val=uh_exc(ne,:)*un;
    elseif x0==x_e(Nelm_e+1)
        unr=basis_1d(Ord_e,-1);
        unl=basis_1d(Ord_e,1);
        Val=(uh_exc(Nelm_e,:)*unl+uh_exc(1,:)*unr)/2;
    end
end

end