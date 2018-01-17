function [ output_args ] = ErrorOrder( err )
    m=size(err,2);
    output_args(1:m-1)=log2(err(1:m-1)./err(2:m));
end

