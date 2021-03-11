%% CFD Final Project: VPE
tic
clear; close all; clc; 
%% Inputs

gamma = 1.4;
R = 287;
T_inf = 300; % Kelvin
P_inf = 101325; % Pa
rho_inf = P_inf/(R*T_inf); % kg/m^3
a_inf = sqrt(gamma*R*T_inf);
M_inf = 0.8; % input
u_inf = M_inf*a_inf;
nsteps = 1000;
omega_vec = [0.5 1 1.5]; % vector of omega values to test

% Part 1 of Mach number correction
M_max = 2.0;
xmcut2 = M_max^2;    %! user specified max Mach number
xmtar2 = xmcut2*(1. + 0.5*(gamma-1.0)*M_inf^2) ...   %! xmach = freestream Mach number
    /(1. + 0.5*(gamma-1.0)*xmcut2);     %! gamma = specific heat ratio

    %   Determines what coordinates to use
    option = 'channel';
    if option == 'channel'
        option = 1;
        file = 'grid.bumps';
    elseif option == 'airfoil'
        option = 2;
        file = 'grid.SD7003';
    end
    
% Read in grid
data = importdata(file);
imax = data(1,1); % finding max i value
jmax = data(1,2); % finding max j value

% Build mesh 
x = reshape(data(2:end,1), [imax,jmax]);
y = reshape(data(2:end,2), [imax,jmax]);


    % Preallocation of wake
    x(1,:) = x(imax-2,:);
    y(1,:) = y(imax-2,:);
    x(imax,:) = x(3,:);
    y(imax,:) = y(3,:);

%% Calculate and store mesh derivatives ( step 3 )
% initial phi equation 
u_inf = M_inf*a_inf; % velocity at infinity 
phi = u_inf.*x;
farfield = phi; % for boundary conditions 

[u,v,ex,ey,zx,zy,xj] = metricDerivsVel(imax,jmax,phi,x,y,M_inf,a_inf,gamma);

%% Flow velocity and Rho Calculation ( step 5 ) (use scrap code for coeffs)
[a,a_spd] = propsCoeffs(u,v,ex,ey,zx,zy,xj,M_inf,imax,jmax);
% calculating initial residual value 

% Residual calculation notes
for j = 2:jmax-1
    for i = 2:imax - 1
        resphi(i,j) = a(i,j,1)*phi(i-1,j-1) + a(i,j,2)*phi(i-1,j) ...
            + a(i,j,3)*phi(i-1,j+1) + a(i,j,4)*phi(i,j-1) ...
            + a(i,j,5)*phi(i,j) + a(i,j,6)*phi(i,j+1) ...
            + a(i,j,7)*phi(i+1,j-1) + a(i,j,8)*phi(i+1,j) ...
            + a(i,j,9)*phi(i+1,j+1);
    end
end

%% Apply Boundary Conditions and recalc flow prop
if option == 1 % channel option
    [phi] = channelBCs(phi,ex,ey,zx,zy,imax,jmax,farfield);
elseif option == 2 % airfoil option
    [phi,x,y] = airfoilBCsmod(phi,ex,ey,zx,zy,imax,jmax,x,y,farfield);
    %   includes wake cut, need to output x and y
end
[u,v,ex,ey,zx,zy,xj] = metricDerivsVel(imax,jmax,phi,x,y,M_inf,a_inf,gamma);
[a,a_spd] = propsCoeffs(u,v,ex,ey,zx,zy,xj,M_inf,imax,jmax);
%% Iterating
for r = 1:length(omega_vec)
    omega = omega_vec(r);
    %omega = 1.2;
    n = 0; % preallocation
    resnorm_rat = 1; % preallocation
    while resnorm_rat > 1e-4 & n <= 5000
        %nsteps = 1000;
        %for n = 1:nsteps
            % calculate residual and sub from freestream 
            n = n+1;
            for j = 2:jmax - 1
                for i = 2: imax - 1
                    res(i,j)= a(i,j,1).*phi(i-1,j-1)+...
                        a(i,j,2).*phi(i-1,j)+...
                        a(i,j,3).*phi(i-1,j+1)+...
                        a(i,j,4).*phi(i,j-1)+...
                        a(i,j,5).*phi(i,j)+...
                        a(i,j,6).*phi(i,j+1)+...
                        a(i,j,7).*phi(i+1,j-1)+...
                        a(i,j,8).*phi(i+1,j)+...
                        a(i,j,9).*phi(i+1,j+1)-...
                        resphi(i,j); 
                end
            end

            % Test for convergence 
            resnorm = 0; 
            for j = 2:jmax - 1
                for i = 2:imax-1
                    resnorm = resnorm + res(i,j)^2; 
                end 
            end
            resnorm = sqrt(abs(resnorm)); 
            if n == 1
                resnorm0 = resnorm; 
            end
            resnorm_rat(n) = resnorm/resnorm0;
            % Storing resnorm_rat's separately
            if omega == omega_vec(1)
                resnorm_rat1(n) = resnorm_rat(n);
            elseif omega == omega_vec(2)
                resnorm_rat2(n) = resnorm_rat(n);
            elseif omega == omega_vec(3)
                resnorm_rat3(n) = resnorm_rat(n);
            end
                
            if resnorm_rat < 1e-4 % checks for convergence 
                break
            end


            % Relaxation Step
            % FSOR 
            for j = 2:jmax - 1
                for i=2:imax-1
                    phi_temp(i,j) = phi(i,j);
                    kernel(i,j) = (resphi(i,j)-(a(i,j,1).*phi(i-1,j-1)+...
                        a(i,j,2).*phi(i-1,j)+...
                        a(i,j,3).*phi(i-1,j+1)+...
                        a(i,j,4).*phi(i,j-1)+...
                        a(i,j,6).*phi(i,j+1)+...
                        a(i,j,7).*phi(i+1,j-1)+...
                        a(i,j,8).*phi(i+1,j)+...
                        a(i,j,9).*phi(i+1,j+1)))./a(i,j,5);
                    phi(i,j) = phi_temp(i,j) + omega*(kernel(i,j)-phi_temp(i,j));
                end
            end

            % BSOR
            for j=jmax-1:-1:2
                for i=imax-1:-1:2
                    phi_temp(i,j) = phi(i,j);
                    kernel(i,j) = (resphi(i,j)-(a(i,j,1).*phi(i-1,j-1)+...
                        a(i,j,2).*phi(i-1,j)+...
                        a(i,j,3).*phi(i-1,j+1)+...
                        a(i,j,4).*phi(i,j-1)+...
                        a(i,j,6).*phi(i,j+1)+...
                        a(i,j,7).*phi(i+1,j-1)+...
                        a(i,j,8).*phi(i+1,j)+...
                        a(i,j,9).*phi(i+1,j+1)))./a(i,j,5);
                    phi(i,j) = phi_temp(i,j) + omega*(kernel(i,j)-phi_temp(i,j));
                end
            end

            % 9d: reapply BCs
            if option == 1
                [phi] = channelBCs(phi,ex,ey,zx,zy,imax,jmax,farfield);
            elseif option == 2
                [phi, x, y] = airfoilBCsmod(phi,ex,ey,zx,zy,imax,jmax,x,y,farfield); % function to apply boundary conditions 
            end
            [u,v,ex,ey,zx,zy,xj] = metricDerivsVel(imax,jmax,phi,x,y,M_inf,a_inf,gamma);
            [a,a_spd] = propsCoeffs(u,v,ex,ey,zx,zy,xj,M_inf,imax,jmax);
            if option == 2
                % Wake Cut Changes 
                phi(1,:) = phi(imax-2,:); 
                u(1,:) = u(imax-2,:);
                v(1,:) = v(imax-2,:);
                %rho(1,:) = rho(imax-2,:);
                phi(imax,:) = phi(3,:);
                u(imax,:) = u(3,:);
                v(imax,:) = v(3,:);
                %rho(imax,:) = rho(3,:);
            end
            resnorm_rat(n) % monitor convergence
            if n >= 30000
                break
            end
    end 
end 
toc
%% Post processing 
tic
v_mag = sqrt(u.^2 + v.^2); % magnitude of velocity 
M_vec = v_mag./a_spd; % mach vector

% Figure 1: Cp vs. x/c
figure(1)
Cp = 1 - (v_mag./(a_inf*M_inf)).^2 ;
plot(x(1:imax/2,1),Cp(1:imax/2,1)) % bottom surface
hold on
plot(x(imax/2:end,1),Cp(imax/2:end,1)) % top surface
hold off
set(gca,'ydir','reverse')% flip y axes
legend('Bottom Surface', 'Top Surface')
xlabel('x/c'); ylabel('Cp'); 
title('Center of Pressure across Channel')
%% Mach over slip surfaces
% Figure 2: Mach over slip surfaces 
figure(2)
plot(x(1:imax/2,1),M_vec(1:imax/2,1));
hold on
plot(x(imax/2:end,1),M_vec(imax/2:end,1));
%{
plot(x(100:200, 1), M_vec(100:200,1), 'k'); hold on; 
plot(flip(x(1:100,1)), flip(M_vec(1:100,1)), 'b'); 
%}
legend('Bottom Surface', 'Top Surface')
xlabel('x/c'); ylabel('M_vec');
title('Mach Along Slip Surface')
%}
%% Figure 3: Contours of velocity magnitude
psi = u.*y - v.*x; % stream function (no rho)
figure(3)
colormap('jet')
contourf(x,y,v_mag, 500,'EdgeColor','none');
hold on
C = caxis; % hold colormap constant
caxis([C(1) C(2)])
contour(x,y,psi,20,'k')
title('Velocity Contours at Omega 1.8')
colorbar; 
%xlim([-1,2]); ylim([-1, 1]); daspect([1 1 1]);
daspect([2 1 1])


%% Figure 4: Convergence
figure(4)
plot(1:length(resnorm_rat1),log10(resnorm_rat1));
hold on
plot(1:length(resnorm_rat2),log10(resnorm_rat2));
plot(1:length(resnorm_rat3),log10(resnorm_rat3));
title('Convergence Plot')
xlabel('Number of Iterations')
ylabel('Convergence Value')
legend('Omega = 0.5','Omega = 1', 'Omega = 1.5');
