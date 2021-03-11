function [a] = coefficients(ex,ey,zx,zy,xj,xme,xmz,rho,imax,jmax)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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

end

