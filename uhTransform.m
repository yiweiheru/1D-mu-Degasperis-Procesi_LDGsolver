function [ uh ] = uhTransform( Nelm,elm_size,U )

uh=zeros(Nelm,elm_size);
for ne=1:Nelm
    for i=1:elm_size
        uh(ne,i)=U((ne-1)*elm_size+i);
    end
end

end

