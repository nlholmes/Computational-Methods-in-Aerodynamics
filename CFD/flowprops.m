function [rho,xmz,xme] = flowprops(u,v,ex,ey,zx,zy,M_inf,imax,jmax)
gamma = 1.4;
R = 287;
T_inf = 300; % Kelvin
P_inf = 101325; % Pa
rho_inf = P_inf/(R*T_inf); % kg/m^3
a_inf = sqrt(gamma*R*T_inf);
u_inf = M_inf*a_inf;

% Calculates flow properties
    for j = 1:jmax
        for i = 1:imax
            % Density initialization
            % Inside the ^(1/(gamma-1)) parenthesis I believe is "term"
            term = 1 + ((gamma - 1)/2)*(M_inf^2 - ((u(i,j)^2 + v(i,j)^2)/a_inf^2));
            %rho(i,j) = rho_inf*(1 + ((gamma - 1)/2)*(M_inf^2 - ((u(i,j)^2 + v(i,j)^2)/a_spd^2)))^(1/(gamma-1));
            %rho(i,j) = rho_inf*(1 + ((gamma - 1)/2)*(M_inf^2 - ((u(i,j)^2 + v(i,j)^2)/a_spd^2)))^(1/(gamma-1));
            rho(i,j) = rho_inf*(term)^(1/(gamma-1));
            % Speed of sound
            % Made a to (i,j)
            a_spd(i,j) = a_inf*sqrt(term);

            % Calculate v dot n
            vdotnxi = (u(i,j)*zx(i,j) + v(i,j)*zy(i,j))/sqrt(zx(i,j)^2 + zy(i,j)^2);
            vdotneta = (u(i,j)*ex(i,j) + v(i,j)*ey(i,j))/sqrt(ex(i,j)^2 + ey(i,j)^2);
            % Calculate local Mach numbers, xmz and xme
            xmz(i,j) = vdotnxi/a_spd(i,j);
            xme(i,j) = vdotneta/a_spd(i,j);

        end
    end
end

