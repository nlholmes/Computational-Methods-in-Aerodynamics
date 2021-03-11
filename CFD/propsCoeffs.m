function [a,a_spd] = propsCoeffs(u,v,ex,ey,zx,zy,xj,M_inf,imax,jmax)
% Calculates flow properties and residual coefficients, half metric derivs
% Flow properties
gamma = 1.4;
R = 287;
T_inf = 300; % Kelvin
P_inf = 101325; % Pa
rho_inf = P_inf/(R*T_inf); % kg/m^3
a_inf = sqrt(gamma*R*T_inf);

%u_inf = M_inf*a_inf;
%{
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


    
    % Initialize for ring of zeros
    a = zeros(imax, jmax, 9);
    % Recalculate rho and coefficients

    for j = 2:jmax-1

        for i = 2:imax-1 % 1 to 2
            rhop(i) = 0.5*(rho(i,j)+rho(i+1,j));
        end
        %rhop(0) = rhop(1);
        rho(1) = rho(2);
        rhop(imax) = rhop(imax-1);

        for i = 1:imax-1 % was 1:imax-1
            xmface = 0.5*(xmz(i,j)+xmz(i+1,j));
            xnui = max(0.0,1.0-1.0/(xmz(i,j)^2 + 1e-8));
            xnuip = max(0.0,1.0-1.0/(xmz(i+1,j)^2 + 1e-8));
            if(xmface > 0.0)
                % added max()
                rhoave = rhop(i) - xnui*(rhop(i)-rhop(max(1,i-1)));
            else
                rhoave = rhop(i) + xnuip*(rhop(i+1)-rhop(i));
            end
            zxave = 0.5*(zx(i,j)+zx(i+1,j));
            zyave = 0.5*(zy(i,j)+zy(i+1,j));
            exave = 0.5*(ex(i,j)+ex(i+1,j));
            eyave = 0.5*(ey(i,j)+ey(i+1,j));
            xjave = 0.5*(xj(i,j)+xj(i+1,j));
            b11 = (zxave^2 + zyave^2)/xjave;
            b12 = (zxave*exave + zyave*eyave)/xjave;
            g11(i) = rhoave*b11;
            g12(i) = rhoave*b12;
        end


        for i = 2:imax-1
            a(i,j,1) =  0.25*g12(i-1);
            a(i,j,3) = -0.25*g12(i-1);
            a(i,j,7) = -0.25*g12(i);
            a(i,j,9) =  0.25*g12(i);
            a(i,j,2) =  g11(i-1);
            a(i,j,8) =  g11(i);
            a(i,j,6) =  0.25*(g12(i) - g12(i-1));
            a(i,j,4) = -0.25*(g12(i) - g12(i-1));
            a(i,j,5) = -(g11(i)+g11(i-1));
        end
    end

    %c ----d/deta(rho*(b12*dphidxi + b22*dphideta))

    for i = 2:imax-1

        for j = 2:jmax-1 % 1 to 2
            rhop(j) = 0.5*(rho(i,j)+rho(i,j+1));
        end
        %rhop(0) = rhop(1);
        rhop(1) = rhop(2);
        rhop(jmax) = rhop(jmax-1);

        for j = 1:jmax-1 % was 1:jmax-1
            xmface = 0.5*(xme(i,j)+xme(i,j+1));
            xnui = max(0.0,1.0-1.0/(xme(i,j)^2 + 1e-8));
            xnuip = max(0.0,1.0-1.0/(xme(i,j+1)^2 + 1e-8));
            if (xmface > 0.0)
                rhoave = rhop(j) - xnui*(rhop(j)-rhop(max(1,j-1)));
            else
                rhoave = rhop(j) + xnuip*(rhop(j+1)-rhop(j));
            end
            zxave = 0.5*(zx(i,j)+zx(i,j+1));
            zyave = 0.5*(zy(i,j)+zy(i,j+1));
            exave = 0.5*(ex(i,j)+ex(i,j+1));
            eyave = 0.5*(ey(i,j)+ey(i,j+1));
            xjave = 0.5*(xj(i,j)+xj(i,j+1));
            b12 = (zxave*exave + zyave*eyave)/xjave;
            b22 = (exave^2 + eyave^2)/xjave;
            h12(j) = rhoave*b12;
            h22(j) = rhoave*b22;
        end

        for j = 2:jmax-1
            a(i,j,1) = a(i,j,1) + 0.25*h12(j-1);
            a(i,j,3) = a(i,j,3) - 0.25*h12(j);
            a(i,j,7) = a(i,j,7) - 0.25*h12(j-1);
            a(i,j,9) = a(i,j,9) + 0.25*h12(j);
            a(i,j,2) = a(i,j,2) - 0.25*(h12(j)-h12(j-1));
            a(i,j,8) = a(i,j,8) + 0.25*(h12(j)-h12(j-1));
            a(i,j,6) = a(i,j,6) + h22(j);
            a(i,j,4) = a(i,j,4) + h22(j-1);
            a(i,j,5) = a(i,j,5) - (h22(j)+h22(j-1));
        end
    end
%}

    % Density calculation
    term = 1 + ((gamma - 1)/2)*(M_inf.^2 - ((u.^2 + v.^2)/a_inf.^2));
    rho = rho_inf*(term).^(1/(gamma-1));
    a_spd = a_inf*sqrt(1+(gamma-1)/2*(M_inf.^2-(u.^2+v.^2)/a_inf.^2)); 
    xmz = ((zx.*u + zy.*v)./sqrt(abs(zx.^2 + zy.^2)))./a_spd; 
    xme = ((ex.*u + ey.*v)./sqrt(abs(ex.^2 + ey.^2)))./a_spd;
    a = zeros(imax, jmax, 9);
    
    % Scrap code for coeffs and metric derivs
    for j = 2:jmax-1
        for i = 1:imax-1 
            rhop(i) = .5*(rho(i,j)+rho(i+1,j));
        end
        rhop(1) = rhop(2);
        rhop(imax) = rhop(imax-1);

        for i = 1:imax-1
            xmface = .5*(xmz(i,j)+xmz(i+1,j));
            xnui = max(0,1-1/(xmz(i,j).^2 + 1e-8));
            xnuip = max(0,1-1/(xmz(i+1,j).^2 + 1e-8));
            if xmface > 0
                rhoave = rhop(i) - xnui*(rhop(i)-rhop(max(1,i-1)));
            else
                rhoave = rhop(i) + xnuip*(rhop(i+1)-rhop(i));
            end
            zxave = .5*(zx(i,j)+zx(i+1,j));
            zyave = .5*(zy(i,j)+zy(i+1,j));
            exave = .5*(ex(i,j)+ex(i+1,j));
            eyave = .5*(ey(i,j)+ey(i+1,j));
            xjave = .5*(xj(i,j)+xj(i+1,j));
            b11 = (zxave.^2 + zyave.^2)/xjave;
            b12 = (zxave*exave + zyave*eyave)/xjave;
            g11(i) = rhoave*b11;
            g12(i) = rhoave*b12;
        end

        for i = 2:imax-1
            a(i,j,1) =  .25*g12(i-1);
            a(i,j,3) = -.25*g12(i-1);
            a(i,j,7) = -.25*g12(i);
            a(i,j,9) =  .25*g12(i);
            a(i,j,2) =  g11(i-1);
            a(i,j,8) =  g11(i);
            a(i,j,6) =  .25*(g12(i) - g12(i-1));
            a(i,j,4) = -.25*(g12(i) - g12(i-1));
            a(i,j,5) = -(g11(i)+g11(i-1));
        end

    end

    % ----d/deta(rho*(b12*dphidxi + b22*dphideta))

    for i = 2:imax-1
        for j = 1:jmax-1 % was 1:jmax-1, no impact noticed
            rhop(j) = .5*(rho(i,j)+rho(i,j+1));
        end
        rhop(1) = rhop(2);
        rhop(jmax) = rhop(jmax-1);

        for j = 1:jmax-1
            xmface = .5*(xme(i,j)+xme(i,j+1));
            xnui = max(0,1-1/(xme(i,j).^2 + 1e-8));
            xnuip = max(0,1-1/(xme(i,j+1).^2 + 1e-8));
            if xmface > 0
                rhoave = rhop(j) - xnui*(rhop(j)-rhop(max(1,j-1)));
            else
                rhoave = rhop(j) + xnuip*(rhop(j+1)-rhop(j));
            end
            zxave = .5*(zx(i,j)+zx(i,j+1));
            zyave = .5*(zy(i,j)+zy(i,j+1));
            exave = .5*(ex(i,j)+ex(i,j+1));
            eyave = .5*(ey(i,j)+ey(i,j+1));
            xjave = .5*(xj(i,j)+xj(i,j+1));
            b12 = (zxave*exave + zyave*eyave)/xjave;
            b22 = (exave.^2 + eyave.^2)/xjave;
            h12(j) = rhoave*b12;
            h22(j) = rhoave*b22;
        end

        for j = 2:jmax-1
            a(i,j,1) = a(i,j,1) + .25*h12(j-1);
            a(i,j,3) = a(i,j,3) - .25*h12(j);
            a(i,j,7) = a(i,j,7) - .25*h12(j-1);
            a(i,j,9) = a(i,j,9) + .25*h12(j);
            a(i,j,2) = a(i,j,2) - .25*(h12(j)-h12(j-1));
            a(i,j,8) = a(i,j,8) + .25*(h12(j)-h12(j-1));
            a(i,j,6) = a(i,j,6) + h22(j);
            a(i,j,4) = a(i,j,4) + h22(j-1);
            a(i,j,5) = a(i,j,5) - (h22(j)+h22(j-1));
        end
    end
    %}
end

