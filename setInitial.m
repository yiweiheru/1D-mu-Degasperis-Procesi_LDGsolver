function [ U0 ] = setInitial(Nelm,elm_size,x,Xexc,uexc)

U0=zeros(Nelm*elm_size,1);
for ne=1:Nelm
    for i=1:elm_size
        if elm_size >= 2
            xtemp = x(ne)+(x(ne+1)-x(ne))*(i-1)/(elm_size-1);
        else
            xtemp = x(ne);
        end
        num=(ne-1)*elm_size+i;
        U0(num,1)=interp1(Xexc,uexc,xtemp);
    end
end

end
