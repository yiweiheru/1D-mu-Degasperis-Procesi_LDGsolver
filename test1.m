% test the condition number of the matrix

% flux_f \in [Dsp,Csv]; flux_q,flux_v \in [C,R,L]
global flux_f flux_q flux_v
flux_f = 'Csv';
flux_q = 'C';
flux_v = 'C';

warning off

Ord = 0;
elm_size = Ord+1;

ir  = 4;
Nelm = 10*2^(ir-1)+1;

period  = 1;
dx = period/Nelm;
x  = 0:dx:period;

[ Amat,Pvmat,Pqmat,massMat,massMat_inv,mu_massMat ] = getAmat(Ord,Nelm,x);

Amat = full(Amat);
Pvmat = full(Pvmat);
Pqmat = full(Pqmat);
massMat = full(massMat);
massMat_inv = full(massMat_inv);
mu_massMat = full(mu_massMat);

format compact

% fprintf("rank of Pqmat*massMat_inv*Pvmat is %i \n", rank(mu_massMat-Amat));
fprintf("rank of Amat is %i \n", rank(Amat));
fprintf("size of Amat is %i \n", size(Amat,1));

fprintf("condition number of Amat is %e \n",  condest(Amat));


% when ord is even & Nelm is odd the condest(Amat) is good.