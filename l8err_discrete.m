function [l8err] = l8err_discrete( Uh,Ue,Ord,Nelm )

elm_size = Ord+1;
uh = uhTransform( Nelm,elm_size,Uh );
ue = uhTransform( Nelm,elm_size,Ue );
u_diff = uh-ue;

NN = 20;
xi = linspace(-1,1,NN+1);
un_at_xi = zeros(elm_size,NN+1);
for k = 1 : NN+1
    un_at_xi(:,k) = basis_1d(Ord,xi(k));
end

l8err = 0;
for ne = 1 : Nelm
    for k = 1 : NN+1
        l8err = max(l8err , abs(u_diff(ne,:)*un_at_xi(:,k)));
    end
end

end
