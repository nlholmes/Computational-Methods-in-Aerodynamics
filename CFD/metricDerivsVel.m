function [u,v,ex,ey,zx,zy,xj] = metricDerivsVel(imax,jmax,phi,x,y,M_inf,a_inf,gamma)
% Inputs:   imax    loop x bound
%           jmax    loop y bound
%           x       mesh x grid
%           y       mesh y grid
% Outputs:  ex
%           ey
%           zx
%           zy
%           xj

% First part of Mach number correction
M_max = 2.0;
xmcut2 = M_max^2;    %! user specified max Mach number
xmtar2 = xmcut2*(1. + 0.5*(gamma-1.0)*M_inf^2) ...   %! xmach = freestream Mach number
    /(1. + 0.5*(gamma-1.0)*xmcut2);     %! gamma = specific heat ratio

    for j = 1:jmax
        for i = 1:imax
            % Eta terms
            if(i == 1)
                ex(i,j) = -(-1.5*y(i,j) + 2*y(i+1,j) - 0.5* y(i+2,j));
                ey(i,j) = (-1.5*x(i,j) + 2*x(i+1,j) - 0.5*x(i+2,j));
                phiz(i,j) = -1.5*phi(i,j) + 2*phi(i+1,j) - 0.5* phi(i+2,j);
            elseif (i == imax)
                ex(i,j) = -(1.5*y(i,j) - 2*y(i-1,j) + 0.5* y(i-2,j));
                ey(i,j) = (1.5*x(i,j) - 2*x(i-1,j) + 0.5*x(i-2,j));
                phiz(i,j) = 1.5*phi(i,j) - 2*phi(i-1,j) + 0.5*phi(i-2,j);
            else
                ex(i,j) = -0.5*(y(i+1,j) - y(i-1,j));
                ey(i,j) = 0.5*(x(i+1,j) - x(i-1,j));
                phiz(i,j) = 0.5*(phi(i+1,j) - phi(i-1,j));
            end

            % Zeta terms
            if(j == 1)
                % moved neg to zy's from zx's
                zx(i,j) = (-1.5*y(i,j) + 2*y(i,j+1) - 0.5* y(i,j+2));
                zy(i,j) = -(-1.5*x(i,j) + 2*x(i,j+1) - 0.5*x(i,j+2));
                phie(i,j) = (-1.5*phi(i,j) + 2*phi(i,j+1) - 0.5*phi(i,j+2));
            elseif (j == jmax)
                zx(i,j) = (1.5*y(i,j) - 2*y(i,j-1) + 0.5* y(i,j-2));
                zy(i,j) = -(1.5*x(i,j) - 2*x(i,j-1) + 0.5*x(i,j-2));
                phie(i,j) = (1.5*phi(i,j) - 2*phi(i,j-1) + 0.5*phi(i,j-2));
            else
                zx(i,j) = 0.5*(y(i,j+1) - y(i,j-1));
                zy(i,j) = -0.5*(x(i,j+1) - x(i,j-1));
                phie(i,j) = 0.5*(phi(i,j+1) - phi(i,j-1));
            end
        end
    end
    xj = zx.*ey - ex.*zy; % calc xj array
    % Velocity calculation
    u = (zx.*phiz + ex.*phie)./xj;
    v = (zy.*phiz + ey.*phie)./xj;
    
    
    % Mach number correction to velocities
    vmagmax = sqrt(xmtar2*a_inf^2);                      %! ainf = freestream sound speed
    vmag = sqrt(u.^2 + v.^2);

    u = min(vmagmax,vmag)./(vmag+1e-6).*u;       %! correction to limit Mach number
    v = min(vmagmax,vmag)./(vmag+1e-6).*v;
end

