function plot_uh( U,Ord,Nelm,x,mark)

uh=uhTransform(Nelm,Ord+1,U);

    x_to_plot = zeros(1,size(x,2));
    uh_to_plot = zeros(1,size(x,2));
    it=0;
    for np=1:size(x,2)
        it=it+1;
        x_to_plot(1,it)=x(np);
        uh_to_plot(1,it)=(evalue_uh(Ord,Nelm,x,uh,x(np)));
    end
if mark == "numerical"
    plot(x_to_plot,uh_to_plot,'^r')
elseif mark =="exact"
    plot(x_to_plot,uh_to_plot,'-b')
end
% if mark == "numerical"
%     x_to_plot = zeros(1,size(x,2));
%     uh_to_plot = zeros(1,size(x,2));
%     it=0;
%     for np=1:size(x,2)
%         it=it+1;
%         x_to_plot(1,it)=x(np);
%         uh_to_plot(1,it)=(evalue_uh(Ord,Nelm,x,uh,x(np)));
%     end
%     plot(x_to_plot,uh_to_plot,'^r')
%     
% elseif mark =="exact"
%     Nelm_plot = 80 ;
%     x_to_plot = zeros(1,Nelm_plot+1);
%     uh_to_plot = zeros(1,Nelm_plot+1);
%     dx_plot = 1/Nelm_plot;
%     x_to_plot = x(1):dx_plot:x(size(x,2));
%     for np = 1 : size(x_to_plot,2)
%         uh_to_plot(np) = evalue_uh(Ord,Nelm,x,uh,x_to_plot(np));
%     end
%     plot(x_to_plot,uh_to_plot,'-b')
% end

end

