function [ex,ey,zx,zy,xj,u,v] = metricDerivsVelocity(imax,jmax,x,y)
% Inputs:   imax    loop x bound
%           jmax    loop y bound
%           x       mesh x grid
%           y       mesh y grid
% Outputs:  ex
%           ey
%           zx
%           zy
%           xj
    for j = 1:jmax
        for i = 1:imax
            % Eta terms
            if(i == 1)
                ex(i,j) = -(-1.5*y(i,j) + 2*y(i+1,j) - 0.5* y(i+2,j));
                ey(i,j) = (-1.5*x(i,j) + 2*x(i+1,j) - 0.5*x(i+2,j));
            elseif (i == imax)
                ex(i,j) = -(1.5*y(i,j) - 2*y(i-1,j) + 0.5* y(i-2,j));
                ey(i,j) = (1.5*x(i,j) - 2*x(i-1,j) + 0.5*x(i-2,j));
            else
                ex(i,j) = -0.5*(y(i+1,j) - y(i-1,j));
                ey(i,j) = 0.5*(x(i+1,j) - x(i-1,j));
            end

            % Zeta terms
            if(j == 1)
                % moved neg to zy's from zx's
                zx(i,j) = (-1.5*y(i,j) + 2*y(i,j+1) - 0.5* y(i,j+2));
                zy(i,j) = -(-1.5*x(i,j) + 2*x(i,j+1) - 0.5*x(i,j+2));
            elseif (j == jmax)
                zx(i,j) = (1.5*y(i,j) - 2*y(i,j-1) + 0.5* y(i,j-2));
                zy(i,j) = -(1.5*x(i,j) - 2*x(i,j-1) + 0.5*x(i,j-2));
            else
                zx(i,j) = 0.5*(y(i,j+1) - y(i,j-1));
                zy(i,j) = -0.5*(x(i,j+1) - x(i,j-1));
            end
        end
    end
    xj = zx.*ey - ex.*zy; % calc xj array
end

