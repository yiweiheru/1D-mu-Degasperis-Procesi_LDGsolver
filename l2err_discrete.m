function [l2err] = l2err_discrete( Uh,Ue,Ord,Nelm )

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

l2err = 0;
for ne = 1 : Nelm
    for k = 1 : NN+1
        l2err = l2err + (u_diff(ne,:)*un_at_xi(:,k))^2;
    end
end
l2err = sqrt(l2err/((NN+1)*Nelm));

end
