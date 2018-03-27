function [ U0 ]  =  setInitial_smo(Nelm,elm_size,x)

U0 = zeros(Nelm*elm_size,1);

for ne = 1:Nelm
    for i = 1:elm_size
        if elm_size >= 2
            xtemp = x(ne)+(x(ne+1)-x(ne))*(i-1)/(elm_size-1);
        else
            xtemp = x(ne);
        end
        num = (ne-1)*elm_size+i;
        % initial condition settings
        U0(num) = 0.1*sin(2*pi*xtemp)+0.5;

        

    end
end
end
