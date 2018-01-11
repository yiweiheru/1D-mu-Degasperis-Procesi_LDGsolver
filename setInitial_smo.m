function [ U0 ]  =  setInitial_smo(Nelm,elm_size,x)

U0 = zeros(Nelm*elm_size,1);

for ne = 1:Nelm
    for i = 1:elm_size

        xtemp = x(ne)+(x(ne+1)-x(ne))*(i-1)/(elm_size-1);
        num = (ne-1)*elm_size+i;
        % initial condition settings
        U0(num) = 0.3*sin(2*pi*xtemp);

        

    end
end
end
