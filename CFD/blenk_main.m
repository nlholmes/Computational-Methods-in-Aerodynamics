%% CFD Final Project Code 
% May 2nd 2019
% Elizabeth Blenk

clear; close all; clc; 

% define known constants
T_i = 300; % Temp initial Kelvin
g = 1.4; % specific density of air
P_i = 101325; % initial pressure 
R = 287; % gas constant of air
rho_i = P_i/(R*T_i); % initial density, I CHANGED THIS 
a_i = sqrt(g*R*T_i); % freestream speed of sound 
M = 0.3; % initial Mach , CHANGE LATER  Mach of 0.18 converges
%converg_p = 1e-4; % given convergence parameter
w = 1.5; % initial value, CHANGE LATER

% Read in the grid file 
data = importdata('grid.SD7003');

imax = data(1,1); % finding max i value
jmax = data(1,2); % finding max j value

% build matrix of x and y points 
x = reshape(data(2:end,1), [imax,jmax]);
y = reshape(data(2:end,2), [imax,jmax]);

% Phi value inially 
phi2 = zeros(imax,jmax); % initialize for speed 
phiF = zeros(imax, jmax); % initialize for speed
phiB = zeros(imax, jmax); % initialize for speed
kernel = zeros(imax, jmax); % initializing for speed
phi_temp = zeros(imax, jmax); % initializing for speed 

% initial phi equation 
Uinf = M*a_i; % velocity at infinity 
for j = 1:jmax
    for i = 1:imax
        phi2(i,j) = Uinf*x(i,j); 
    end
end
farfield = phi2; % for boundary conditions 

% calling funciton finds all the derivative things 
[ex, ey, zx, zy, xj2, u2, v2, rho2, a, speed_s] = met_der(x, y, imax, jmax, phi2, M, a_i, g, rho_i);

% calculating initial residual value 
res_0 = zeros(imax, jmax); % preallocating for speed
res = res_0; % preallocating for speed 

for j = 2:jmax - 1
    for i = 2: imax - 1
        res0(i,j) = a(i,j,1).*phi2(i-1,j-1)+...
            a(i,j,2).*phi2(i-1,j)+...
            a(i,j,3).*phi2(i-1,j+1)+...
            a(i,j,4).*phi2(i,j-1)+...
            a(i,j,5).*phi2(i,j)+...
            a(i,j,6).*phi2(i,j+1)+...
            a(i,j,7).*phi2(i+1,j-1)+...
            a(i,j,8).*phi2(i+1,j)+...
            a(i,j,9).*phi2(i+1,j+1);
    end
end
resphi2 = res0;

%% Apply Boundary Conditions and recalc flow prop
[phi2, x, y] = boundary_airfoil(x, y, phi2, ex, ey, zx, zy, imax, jmax, farfield); % function to apply boundary conditions 
[ex, ey, zx, zy, xj, u, v, rho, a, speed_s] = met_der(x, y, imax, jmax, phi2, M, a_i, g, rho_i); % find derivatives

%% Now to run through looping things 
converge = 0;  
n = 0; % counter variable 
%while converge == 0 
for h = 1:100
    n = n + 1; % number of iterations 
    
    % calculate residual and sub from freestream 
    for j = 2:jmax - 1
        for i = 2: imax - 1
            res(i,j)= a(i,j,1).*phi2(i-1,j-1)+...
                a(i,j,2).*phi2(i-1,j)+...
                a(i,j,3).*phi2(i-1,j+1)+...
                a(i,j,4).*phi2(i,j-1)+...
                a(i,j,5).*phi2(i,j)+...
                a(i,j,6).*phi2(i,j+1)+...
                a(i,j,7).*phi2(i+1,j-1)+...
                a(i,j,8).*phi2(i+1,j)+...
                a(i,j,9).*phi2(i+1,j+1)-...
                res0(i,j); 
        end
    end
    
    % Test for convergence 
    resnorm = 0; 
    for j = 2:jmax - 1
        for i = 2:imax-1
            resnorm = resnorm + res(i,j).*res(i,j); 
        end 
    end
    resnorm = sqrt(abs(resnorm)); 
    if n == 1
        resnorm0 = resnorm; 
    end
    %if resnorm/resnorm0 < converg_p % checks for convergence 
        %break
    %end
    %p_val(n) = resnorm/resnorm0;  % plotting the convergence history 
    %p2_val(n) = n; 
    
    % Relaxation Step
    % FSOR 
    for j = 2:jmax - 1
        for i=2:imax-1
            phi_temp(i,j) = phi2(i,j);
            kernel(i,j) = (res0(i,j)-(a(i,j,1).*phi2(i-1,j-1)+...
                a(i,j,2).*phi2(i-1,j)+...
                a(i,j,3).*phi2(i-1,j+1)+...
                a(i,j,4).*phi2(i,j-1)+...
                a(i,j,6).*phi2(i,j+1)+...
                a(i,j,7).*phi2(i+1,j-1)+...
                a(i,j,8).*phi2(i+1,j)+...
                a(i,j,9).*phi2(i+1,j+1)))./a(i,j,5);
            phi2(i,j) = phi_temp(i,j)+w*(kernel(i,j)-phi_temp(i,j));
        end
    end
    
    % BSOR
    for j=jmax-1:-1:2
        for i=imax-1:-1:2
            phi_temp(i,j) = phi2(i,j);
            kernel(i,j) = (res0(i,j)-(a(i,j,1).*phi2(i-1,j-1)+...
                a(i,j,2).*phi2(i-1,j)+...
                a(i,j,3).*phi2(i-1,j+1)+...
                a(i,j,4).*phi2(i,j-1)+...
                a(i,j,6).*phi2(i,j+1)+...
                a(i,j,7).*phi2(i+1,j-1)+...
                a(i,j,8).*phi2(i+1,j)+...
                a(i,j,9).*phi2(i+1,j+1)))./a(i,j,5);
            phi2(i,j) = phi_temp(i,j)+w*(kernel(i,j)-phi_temp(i,j));
        end
    end
    
    % Apply Boundary Conditions finding new values for things 
    [phi2, x, y] = boundary_airfoil(x, y, phi2, ex, ey, zx, zy, imax, jmax, farfield); % function to apply boundary conditions 
    [ex, ey, zx, zy, xj, u, v, rho, a, speed_s] = met_der(x, y, imax, jmax, phi2, M, a_i, g, rho_i);
    
    % Wake Cut Changes 
    phi2(1,:) = phi2(imax-2,:); 
    u(1,:) = u(imax-2,:);
    v(1,:) = v(imax-2,:);
    rho(1,:) = rho(imax-2,:);
    phi2(imax,:) = phi2(3,:);
    u(imax,:) = u(3,:);
    v(imax,:) = v(3,:);
    rho(imax,:) = rho(3,:);
end 
%% Post processing 
v_m = sqrt(u.^2 + v.^2); % magnitude of velocity 
%M_vec = v_m/real(speed_s); % mach vector

% Figure 1: Cp vs. x/c along slip surfaces (just for airfoil) 
%{
figure
Cp = 1 - (v./(a_i*M)).^2 ; 
plot(x(100:200),Cp(100:200)); hold on; % these plot the x and y values of the airfoil 
plot(flip(x(1:100,1)), flip(Cp(1:100,1))); 
legend('Bottom Surface', 'Top Surface','Location', 'SouthEast')
xlabel('x'); ylabel('Cp'); 
title('Center of Pressure along Airfoil')
%}
% Figure 2: Mach # along slip surfaces 
%{
figure 
plot(x(100:200, 1), M_vec(100:200,1), 'k'); hold on; 
plot(flip(x(1:100,1)), flip(M_vec(1:100,1)), 'b'); 
legend('Bottom Surface', 'Top Surface', 'Location', 'SouthEast')
xlabel('x'); ylabel('y');
title('Mach Along Slip Surface')
%}
% Figure 3: Contours of veloity magnitude
figure 
contourf(x(:,10:end-10),y(:,10:end-10),v_m(:,10:end-10), 500,'EdgeColor','none'); colorbar; colormap('jet')
xlim([-1,1.5]); ylim([-1, 1]); daspect([1 1 1]);

% Figure 4: 
%figure
%plot(log10(p2_val), p_val); hold on
%title('Convergence Plot')
%xlabel('Number of Iterations')
%ylabel('Convergence Value')

% BOUNDARY CONDITIONS FOR THIS ONE SEEM TOO HIGH FOR SOME REASON TRY TAKING
% SOME VALUES OUT OF PLOT??
%%
phi2 = phi2;
u2 = u;
v2 = v;
rho2 = rho;