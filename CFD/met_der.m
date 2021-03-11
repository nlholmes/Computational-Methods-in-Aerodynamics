function [ex, ey, zx, zy, xj, u, v, rho, a, speed_s] = met_der(x, y, imax, jmax, phi, M, a_i, g, rho_i)
% subroutine to calculate metric derivatives and inverse Jacobians for a X-Y mesh. 
%input should include:
%    x(i,j) and y(i,j) arrays that store the mesh coordinates

% ADDED IN: Mach correction part 1
% No real change noticed, but may aid in higher Mach calculations
% Actually, seems to have removed complex values from calculations
M_max = 2.0;
xmcut2 = M_max^2;    %! user specified max Mach number
xmtar2 = xmcut2*(1. + 0.5*(g-1.0)*M^2) ...   %! xmach = freestream Mach number
    /(1. + 0.5*(g-1.0)*xmcut2);     %! gamma = specific heat ratio

% initializing values for speed
 
a = zeros(imax, jmax, 9);   
ex = zeros(imax, jmax); 
ey = zeros(imax, jmax); 
zx = zeros(imax, jmax); 
zy = zeros(imax, jmax);  
phi_e = zeros(imax, jmax); 
phi_z = zeros(imax, jmax); 
u = zeros(imax, jmax); 
v = zeros(imax, jmax); 


% Loop for finding metric derivatives 
for j = 1 : jmax
    for i = 1 : imax
        % defining eta for x and y
        if i == 1
            ex(i,j) = -(-1.5*y(i,j) + 2*y(i+1,j) - .5*y(i+2,j));
            ey(i,j) = (-1.5*x(i,j) + 2*x(i+1,j) - .5*x(i+2,j));
            phi_z(i,j) = (-1.5*phi(i,j)+2*phi(i+1,j)-.5*phi(i+2,j));
        elseif i == imax
            ex(i,j) = -(1.5*y(i,j) - 2*y(i-1,j) + .5*y(i-2,j));
            ey(i,j) = (1.5*x(i,j) - 2*x(i-1,j) + .5*x(i-2,j));
            phi_z(i,j)=(1.5*phi(i,j)-2*phi(i-1,j)+.5*phi(i-2,j));
        else
            ex(i,j) = -0.5*(y(i+1,j) - y(i-1,j));
            ey(i,j) = 0.5*(x(i+1,j) - x(i-1,j));
            phi_z(i,j)= 0.5*(phi(i+1,j)-phi(i-1,j));
        end
        
        % defining xi for x and y
        if j == 1
            zx(i,j) = (-1.5*y(i,j) + 2*y(i,j+1) - .5*y(i,j+2));
            zy(i,j) = -(-1.5*x(i,j) + 2*x(i,j+1) - .5*x(i,j+2));
            phi_e(i,j) = (-1.5*phi(i,j)+2*phi(i,j+1)-.5*phi(i,j+2));
        elseif j == jmax
            zx(i,j) = (1.5*y(i,j) - 2*y(i,j-1) + .5*y(i,j-2)); % check these signs 
            zy(i,j) = -(1.5*x(i,j) - 2*x(i,j-1) + .5*x(i,j-2)); % check these signs 
            phi_e(i,j) = (1.5*phi(i,j)-2*phi(i,j-1)+.5*phi(i,j-2)); % check these signs 
        else
            zx(i,j) = 0.5*(y(i,j+1) - y(i,j-1));
            zy(i,j) = -0.5*(x(i,j+1) - x(i,j-1));
            phi_e(i,j) = 0.5*(phi(i,j+1)-phi(i,j-1));
        end
        
    end
end
xj = zx.*ey - ex.*zy;
u = (zx.*phi_z+ex.*phi_e)./xj;
v = (zy.*phi_z+ey.*phi_e)./xj;


    % Mach number correction part 2(ADDED IN)
    vmagmax = sqrt(xmtar2*a_i^2);                      %! ainf = freestream sound speed
    vmag = sqrt(u.^2 + v.^2);
    u = min(vmagmax,vmag)./(vmag+1e-6).*u;       %! correction to limit Mach number
    v = min(vmagmax,vmag)./(vmag+1e-6).*v;

    
% finding rho 
rho = rho_i*(1 + ((g - 1)/2)*(M.^2 - ((u.^2 + v.^2)/a_i.^2))).^(1/(g-1));
speed_s = a_i*sqrt(1+(g-1)/2*(M.^2-(u.^2+v.^2)/a_i.^2)); 
xmz = ((zx.*u + zy.*v)./sqrt(abs(zx.^2 + zy.^2)))./speed_s; 
xme = ((ex.*u + ey.*v)./sqrt(abs(ex.^2 + ey.^2)))./speed_s;

% Jack's scrap code for calculating metric derivatives 
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
    for j = 2:jmax-1 % was 1:jmax-1, no impact noticed
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

end