function [ U_final ] = get_exac( elm_size,Nelm,x,Ord_e,ir_e,Uexc,period )

Nelm_e = 10*2^(ir_e-1);
elm_size_e = Ord_e+1;
x_e = 0:(period/Nelm_e):period;
uh_exc = uhTransform( Nelm_e,elm_size_e,Uexc );
U_final=zeros(Nelm*elm_size,1);
for ne=1:Nelm
    for i=1:elm_size
        xtemp=x(ne)+(x(ne+1)-x(ne))*(i-1)/(elm_size-1);
        num=(ne-1)*elm_size+i;
        U_final(num,1) = evalue_uh( Ord_e,Nelm_e,x_e,uh_exc,xtemp );
    end
end


end

