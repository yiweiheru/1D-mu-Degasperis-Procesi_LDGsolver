function [ NFmat ] = getAmat_flux(Nelm,elm_size,Fmat,flux_type)
% This func. give different fluxes at nodes.

nei=zeros(Nelm,2);
for ne=1:Nelm
    if ne==1
        nei(ne,1)=Nelm;
        nei(ne,2)=2;
    elseif ne==Nelm
        nei(ne,1)=Nelm-1;
        nei(ne,2)=1;
    else
        nei(ne,1)=ne-1;
        nei(ne,2)=ne+1;
    end
end

switch flux_type
    case 'L'
        % In this case: flux are chosen as upwind u^-
        isp  = zeros(2*Nelm * elm_size^2 + 1,1);
        jsp  = zeros(2*Nelm * elm_size^2 + 1,1);
        sdat = zeros(2*Nelm * elm_size^2 + 1,1);
        it = 0;
        for ne = 1:Nelm
            for i = 1:elm_size
                for j = 1:elm_size
                    it = it+1;
                    isp(it)  = (ne-1)*elm_size+i;
                    jsp(it)  = (nei(ne,1)-1)*elm_size+j;
                    sdat(it) = Fmat(i,j,1,2,ne);
                    
                    it=it+1;
                    isp(it)  = (ne-1)*elm_size+i;
                    jsp(it)  = (ne-1)*elm_size+j;
                    sdat(it) = -Fmat(i,j,2,2,ne);
                end
            end
        end
        it = it+1;
        isp(it)  = Nelm*elm_size;
        jsp(it)  = Nelm*elm_size;
        sdat(it) = 0;
        NFmat = sparse(isp,jsp,sdat);
        
    case 'R'
        % In this case: flux are chosen as downwind u^+
        isp  = zeros(2*Nelm * elm_size^2 + 1,1);
        jsp  = zeros(2*Nelm * elm_size^2 + 1,1);
        sdat = zeros(2*Nelm * elm_size^2 + 1,1);
        it = 0;
        for ne = 1:Nelm
            for i = 1:elm_size
                for j = 1:elm_size
                    it = it+1;
                    isp(it)  = (ne-1)*elm_size+i;
                    jsp(it)  = (ne-1)*elm_size+j;
                    sdat(it) = Fmat(i,j,1,1,ne);
                    
                    it=it+1;
                    isp(it)  = (ne-1)*elm_size+i;
                    jsp(it)  = (nei(ne,2)-1)*elm_size+j;
                    sdat(it) = -Fmat(i,j,2,1,ne);
                end
            end
        end
        it = it+1;
        isp(it)  = Nelm*elm_size;
        jsp(it)  = Nelm*elm_size;
        sdat(it) = 0;
        NFmat = sparse(isp,jsp,sdat);
        
    case 'C'
        % In this case: flux are chosen as central 0.5*(u^+ + u^-)
        isp  = zeros(2*2*Nelm * elm_size^2 + 1,1);
        jsp  = zeros(2*2*Nelm * elm_size^2 + 1,1);
        sdat = zeros(2*2*Nelm * elm_size^2 + 1,1);
        it = 0;
        for ne = 1:Nelm
            for i = 1:elm_size
                for j = 1:elm_size
                    it = it+1;
                    isp(it)  = (ne-1)*elm_size+i;
                    jsp(it)  = (ne-1)*elm_size+j;
                    sdat(it) = 0.5*Fmat(i,j,1,1,ne);
           
                    it = it+1;
                    isp(it)  = (ne-1)*elm_size+i;
                    jsp(it)  = (nei(ne,1)-1)*elm_size+j;
                    sdat(it) = 0.5*Fmat(i,j,1,2,ne);
                    
                    it = it+1;
                    isp(it)  = (ne-1)*elm_size+i;
                    jsp(it)  = (nei(ne,2)-1)*elm_size+j;
                    sdat(it) = -0.5*Fmat(i,j,2,1,ne);
                    
                    it = it+1;
                    isp(it)  = (ne-1)*elm_size+i;
                    jsp(it)  = (ne-1)*elm_size+j;
                    sdat(it) = -0.5*Fmat(i,j,2,2,ne);
                end
            end
        end
        it = it+1;
        isp(it)  = Nelm*elm_size;
        jsp(it)  = Nelm*elm_size;
        sdat(it) = 0;
        NFmat = sparse(isp,jsp,sdat);
        
end

