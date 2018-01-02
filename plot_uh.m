function [x_to_plot,uh_to_plot]=plot_uh( U,Ord,Nelm,x)

uh=uhTransform(Nelm,Ord+1,U);

x_to_plot = zeros(1,size(x,2));
uh_to_plot = zeros(1,size(x,2));
it=0;
for np=1:size(x,2)
    it=it+1;
    x_to_plot(1,it)=x(np);
    uh_to_plot(1,it)=(evalue_uh(Ord,Nelm,x,uh,x(np)));
end

end

