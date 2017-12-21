function [ Val ] = evalue_uh( Ord,Nelm,x,uh,x0 )

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
    if x0==x(ne) 
        unr=basis_1d(Ord,-1);
        unl=basis_1d(Ord,1);
        Val=(uh(ne,:)*unr+uh(nei(ne,1),:)*unl)/2;
    elseif x0>x(ne) && x0<x(ne+1)
        xi=(x0-x(ne))/(x(ne+1)-x(ne))*2+(-1);
        un=basis_1d(Ord,xi);
        Val=uh(ne,:)*un;
    elseif x0==x(Nelm+1)
        unr=basis_1d(Ord,-1);
        unl=basis_1d(Ord,1);
        Val=(uh(Nelm,:)*unl+uh(1,:)*unr)/2;
    end
end

end