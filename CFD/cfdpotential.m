function [u,v] = cfdpotential(ex,ey,zx,zy,xj,imax,jmax,phi,xmtar2,a_inf)
% Inputs:   imax    loop x bound
%           jmax    loop y bound
%           phi     matrix of phi values
% Outputs:  phiz    flow potential in computational coordinate zeta direction
%           phie    flow potential in computational coordinate eta direction


    for j = 1:jmax
        for i = 1:imax
            % Zeta terms
            if(i == 1)
                phiz(i,j) = -1.5*phi(i,j) + 2*phi(i+1,j) - 0.5* phi(i+2,j);
            elseif (i == imax)
                phiz(i,j) = 1.5*phi(i,j) - 2*phi(i-1,j) + 0.5*phi(i-2,j);
            else
                phiz(i,j) = 0.5*(phi(i+1,j) - phi(i-1,j));
            end

            % Eta terms
            % Removed negatives from the front of all 3 
            if(j == 1)
                phie(i,j) = (-1.5*phi(i,j) + 2*phi(i,j+1) - 0.5*phi(i,j+2));
            elseif (j == jmax)
                phie(i,j) = (1.5*phi(i,j) - 2*phi(i,j-1) + 0.5*phi(i,j-2));
            else
                phie(i,j) = 0.5*(phi(i,j+1) - phi(i,j-1));
            end
        end
    end
    
    u = (zx.*phiz + ex.*phie)./xj;
    v = (zy.*phiz + ey.*phie)./xj;
    
    
    % Mach number correction
    vmagmax = sqrt(xmtar2*a_inf^2);                      %! ainf = freestream sound speed
    vmag = sqrt(u.^2 + v.^2);

    u = min(vmagmax,vmag)./(vmag+1e-6).*u;       %! correction to limit Mach number
    v = min(vmagmax,vmag)./(vmag+1e-6).*v;
    
end

