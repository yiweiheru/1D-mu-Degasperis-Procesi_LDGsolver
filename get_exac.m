function [ U_final ] = get_exac( elm_size,Nelm,x,Ord_e,ir_e,x_e,Nelm_e,elm_size_e,Uexc )

uh_exc = uhTransform( Nelm_e,elm_size_e,Uexc );
U_final= zeros(Nelm*elm_size,1);
for ne=1:Nelm
    for i=1:elm_size
        if elm_size >= 2
            xtemp = x(ne)+(x(ne+1)-x(ne))*(i-1)/(elm_size-1);
        else
            xtemp = x(ne);
        end
        num=(ne-1)*elm_size+i;
        U_final(num,1) = evalue_uh( Ord_e,Nelm_e,x_e,uh_exc,xtemp );
    end
end


end

