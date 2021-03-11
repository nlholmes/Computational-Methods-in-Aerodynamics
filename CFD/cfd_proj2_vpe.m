%% CFD Project 2: Solution of VPE in Generalized Coordinates, airfoil in tunnel
clear all;close all;clc;
%% Todo
% 1.    Read in input parameters
% 2.    Read in mesh
% 3.    Calculate and store mesh derivatives and cell area
% 4.    Initialize the flow
% 5.    Calculate u,v,rho,a
% 6.    Calculate initial residual
% 7.    Apply boundary conditions
% 8.    Recalculate flow properties and a again (step 5)
% 9.    Iterate
% 10.   Post-process

%% Inputs
%u_inf = 200;
%v_inf = 200;
gamma = 1.4;
R = 287;
T_inf = 300; % Kelvin
P_inf = 101325; % Pa
rho_inf = P_inf/(R*T_inf); % kg/m^3
a_inf = sqrt(gamma*R*T_inf);
M_inf = 0.2; % input
u_inf = M_inf*a_inf;
nsteps = 1000;
omega = 1.8;

% Part of Mach number correction
M_max = 2.0;
xmcut2 = M_max^2;    %! user specified max Mach number
xmtar2 = xmcut2*(1. + 0.5*(gamma-1.0)*M_inf^2) ...   %! xmach = freestream Mach number
    /(1. + 0.5*(gamma-1.0)*xmcut2);     %! gamma = specific heat ratio

% To make the code adaptive
%option = 'channel'; % 'airfoil' or 'channel
option = 'airfoil';
if option == 'channel'
    option = 1;
elseif option == 'airfoil'
    option = 2;
end
%% Read in input and mesh ( step 1 and step 2 )
channel = importdata('grid.bumps');
airfoil = importdata('grid.SD7003');
chanimax = channel(1,1);
chanjmax = channel(1,2);
foilimax = airfoil(1,1);
foiljmax = airfoil(1,2);

% Reshape data into meshgrids
% Channel
channelX = reshape(channel(2:end,1), [201,51]);
channelY = reshape(channel(2:end,2), [201,51]);
% Airfoil
airfoilX = reshape(airfoil(2:end,1), [200,81]);
airfoilY = reshape(airfoil(2:end,2), [200,81]);

%plot(channelX,channelY)
%plot(airfoilX,airfoilY)
%% Calculate and store mesh derivatives ( step 3 )
% for now deal with airfoil
% Channel turns out like u and v plots
% Something is definitely wrong with the code.
if option == 1
    imax = chanimax;
    jmax = chanjmax;
    x = channelX;
    y = channelY;
    phi = u_inf.*x;
elseif option == 2
    imax = foilimax;
    jmax = foiljmax;
    x = airfoilX;
    y = airfoilY;
    phi = u_inf.*x;
    
    % Initialize Wake Cut
    %{ 
    % without initializing and in boundary, stronger wake
    x(1,:) = x(imax-2,:);
    y(1,:) = y(imax-2,:);
    x(imax,:) = x(3,:);
    y(imax,:) = y(3,:);
    %}
    %{
    phi(1,:) = phi(imax-2,:);
    
    u(1,:) = u(imax-2,:);
    v(1,:) = v(imax-2,:);
    rho(1,:) = rho(imax-2,:);
    phi(imax,:) = phi(3,:);
    
    u(imax,:) = u(3,:);
    v(imax,:) = v(3,:);
    rho(imax,:) = rho(3,:);
    %}
end



% Previous code (HW 4) put into metric derivative function, altered
[ex,ey,zx,zy,xj] = metricDerivs(imax,jmax,x,y);
%% Flow initialization ( step 4 )
farfield = phi; % farfield BC initializer, later use in the BC functions

% Density initialization
% use cfdpotential user function to get phiz and phie for velocity component calculation
% Calculate velocity components in grid
[u,v] = cfdpotential(ex,ey,zx,zy,xj,imax,jmax,phi,xmtar2,a_inf);

[rho,xmz,xme] = flowprops(u,v,ex,ey,zx,zy,M_inf,imax,jmax);


%% Flow velocity and Rho Calculation ( step 5 ) (use scrap code for coeffs)

a = coefficients(ex,ey,zx,zy,xj,xme,xmz,rho,imax,jmax);

%% Calculate Initial Residual ( step 6 )
% Residual calculation notes
% do in double loop (2:max-1) for both i and j
% res = (phi(i+1,j) - 2*phi(i,j)@k+1 + phi(i-1,j))/dx^2 + (phi(i,j+1) -
% 2*phi(i,j)@k+1 + phi(i,j-1))/dy^2 - f(i,j) tends to zero
% xL2norm = xL2norm + res^2
% can plot log10(xL2rat vs k)
% stop loop once res is 0.0001 or less
%   if log10(xL2rat) < 10e-4, stop loop
% This was right before classes of relaxation methods, talking about pt
% jacobi and etc

% From later notes on the flow chart: step 9a
% In notes chris sent, the residual calc did not have the initial residual
% (resphi(i,j)), so dont include it here, but do so in step 9
for j = 2:jmax-1
    for i = 2:imax - 1
        resphi(i,j) = a(i,j,1)*phi(i-1,j-1) + a(i,j,2)*phi(i-1,j) ...
            + a(i,j,3)*phi(i-1,j+1) + a(i,j,4)*phi(i,j-1) ...
            + a(i,j,5)*phi(i,j) + a(i,j,6)*phi(i,j+1) ...
            + a(i,j,7)*phi(i+1,j-1) + a(i,j,8)*phi(i+1,j) ...
            + a(i,j,9)*phi(i+1,j+1);
    end
end
resphi1 = resphi;
%}

%% Apply Boundary Conditions ( step 7 )
% Recent BC notes
% Tangency, Far field, wake cut

if option == 1 % channel option
    [phi] = channelBCs(phi,ex,ey,zx,zy,imax,jmax,farfield);
elseif option == 2 % airfoil option
    [phi,x,y] = airfoilBCs(phi,ex,ey,zx,zy,imax,jmax,x,y,farfield);
    %   includes wake cut, need to output x and y
end


%% Recalculate Flow Properties ( step 8 )
% Recent notes/repeat code from step 5

% Recalculate new u and v components from new phi values
[u,v] = cfdpotential(ex,ey,zx,zy,xj,imax,jmax,phi,xmtar2,a_inf);

% Recalculate vdotn (flow props)
[rho,xmz,xme] = flowprops(u,v,ex,ey,zx,zy,M_inf,imax,jmax);

% Recalculate rho and coefficients
[a] = coefficients(ex,ey,zx,zy,xj,xme,xmz,rho,imax,jmax);

%% Iterate ( step 9 )
for n = 1:nsteps
        resnorm = 0;
        for j = 2:jmax-1
            for i = 2:imax - 1
                % Recalculate residual
                res(i,j) = a(i,j,1)*phi(i-1,j-1) + a(i,j,2)*phi(i-1,j) ...
                    + a(i,j,3)*phi(i-1,j+1) + a(i,j,4)*phi(i,j-1) ...
                    + a(i,j,5)*phi(i,j) + a(i,j,6)*phi(i,j+1) ...
                    + a(i,j,7)*phi(i+1,j-1) + a(i,j,8)*phi(i+1,j) ...
                    + a(i,j,9)*phi(i+1,j+1) - resphi(i,j); % assume resphi is initial residual
                resnorm = resnorm + res(i,j)^2; % residual norm
            end
        end
        resnorm = sqrt(resnorm); % squareroot of residual
        if n == 1
            resnormphi = resnorm; % resnormphi is the first step resnorm calc
        end
        
        if (resnorm/resnormphi < 1e-4) % resphi is initial residual norm and initial residual
            break
        end

    % Begin SSOR (FSOR + BSOR)
    % FSOR
    for j = 2:jmax-1 % 2
        for i = 2:imax-1 % 2
            phitmp = phi(i,j);
            kernel = (resphi(i,j) - (a(i,j,1)*phi(i-1,j-1) + a(i,j,2)*phi(i-1,j) ...
                + a(i,j,3)*phi(i-1,j+1) + a(i,j,4)*phi(i,j-1) ...
                + a(i,j,6)*phi(i,j+1) ...
                + a(i,j,7)*phi(i+1,j-1) + a(i,j,8)*phi(i+1,j) ...
                + a(i,j,9)*phi(i+1,j+1)))/a(i,j,5); 
            phi(i,j) = phitmp + omega*(kernel - phitmp);
        end
    end
    % BSOR: is the same body but the i-j loop order traverse reversed
    
    for j = jmax-1:-1:2
        for i = imax-1:-1:2
            phitmp = phi(i,j);
            % missing parenth
            kernel = (resphi(i,j) - (a(i,j,1)*phi(i-1,j-1) + a(i,j,2)*phi(i-1,j) ...
                + a(i,j,3)*phi(i-1,j+1) + a(i,j,4)*phi(i,j-1) ...
                + a(i,j,6)*phi(i,j+1) ...
                + a(i,j,7)*phi(i+1,j-1) + a(i,j,8)*phi(i+1,j) ...
                + a(i,j,9)*phi(i+1,j+1)))/a(i,j,5); 
            phi(i,j) = phitmp + omega*(kernel - phitmp);
        end
    end
    
    % 9d: reapply BCs
    if option == 1
        [phi] = channelBCs(phi,ex,ey,zx,zy,imax,jmax,farfield);
    elseif option == 2
        [phi,x,y] = airfoilBCs(phi,ex,ey,zx,zy,imax,jmax,x,y,farfield);
    end
    
    %9e: recalculate flow properties u,v,rho,a_spd(i,j) (maybe coefficients)
    [u,v] = cfdpotential(ex,ey,zx,zy,xj,imax,jmax,phi,xmtar2,a_inf);
    
    [rho,xmz,xme] = flowprops(u,v,ex,ey,zx,zy,M_inf,imax,jmax);
    

    % Redoing coefficient calculations
    [a] = coefficients(ex,ey,zx,zy,xj,xme,xmz,rho,imax,jmax);
    
    % Apply wake cut
    
    if option == 2
        phi(1,:) = phi(imax-2,:);
        u(1,:) = u(imax-2,:);
        v(1,:) = v(imax-2,:);
        rho(1,:) = rho(imax-2,:);
        phi(imax,:) = phi(3,:);
        u(imax,:) = u(3,:);
        v(imax,:) = v(3,:);
        rho(imax,:) = rho(3,:);
    end
    
    
    resnorm/resnormphi % convergence notifier
    
end

%% Post-processing
% Recall: in HW 4 to plot u and v, had to subtract initial value to get the
% minute differences across the field

vel_mag = sqrt(u.^2 + v.^2);
% Subtract initial value of u_inf
figure(2)
colormap('jet')
% had u-u_inf here, with the small number added at the division by zero of
% a(i,j,5), the plot looks...good?
localMach = sqrt(xmz.^2 + xme.^2);
%contourf(x,y,localMach,100,'edgecolor','none');
% first three rows are wake cut
%contourf(x,y,vel_mag,500,'edgecolor','none');
contourf(x,y,vel_mag,500,'edgecolor','none')
if option == 1
    daspect([2 1 1])
elseif option == 2
    axis([-0.5 1.5 -0.5 0.5])
end
colorbar;
%% Coefficient scrap code, and fixed (from email)
%{
c **************************************************************
c ----- scrapcode for calculating matrix coefficients
c **************************************************************

% NEED
% What is xmz?
%   local mach in zeta direction
%   xme is local mach in eta direction
% Define/initialize rho

for j = 2:jmax-1

    for i = 1:imax-1
        rhop(i) = 0.5*(rho(i,j)+rho(i+1,j));
    end
    rhop(0) = rhop(1);
    rhop(imax) = rhop(imax-1);

    for i = 1:imax-1
        xmface = 0.5*(xmz(i,j)+xmz(i+1,j));
        xnui = max(0.0,1.0-1.0/(xmz(i,j)^2 + 1e-8));
        xnuip = max(0.0,1.0-1.0/(xmz(i+1,j)^2 + 1e-8));
        if(xmface > 0.0)
            rhoave = rhop(i) - xnui*(rhop(i)-rhop(i-1));
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
%end

%c ----d/deta(rho*(b12*dphidxi + b22*dphideta))

for i = 2:imax-1

    for j = 1:jmax-1
        rhop(j) = 0.5*(rho(i,j)+rho(i,j+1));S
        enddo
        rhop(0) = rhop(1);
        rhop(jmax) = rhop(jmax-1);

        for j = 1:jmax-1
            xmface = 0.5*(xme(i,j)+xme(i,j+1));
            xnui = max(0.0,1.0-1.0/(xme(i,j)^2 + 1e-8));
            xnuip = max(0.0,1.0-1.0/(xme(i,j+1)^2 + 1e-8));
            if (xmface > 0.0)
                rhoave = rhop(j) - xnui*(rhop(j)-rhop(j-1));
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

% return
%end
%}
%% More scrapcode
%{

c **************************************************************
c ----- scrapcode for applying surface tangency bcs
c ----- set V-dot-n = 0 =>b12*dphidxi + b22*dphideta = 0.0
c **************************************************************

c ---- surface at j=1

      j=1
      do i=2,ii-1
       dphidxi(i) = 0.5*(phi(i+1,j)-phi(i-1,j))  !lagged derivative
      enddo

      do i=2,ii-1
       b12 = (zx(i,j)*ex(i,j) + zy(i,j)*ey(i,j))
       b22 = (ex(i,j)**2 + ey(i,j)**2)
       phi(i,j) = (-b12*dphidxi(i) - (2.0*phi(i,j+1)-0.5*phi(i,j+2))*b22)/(-1.5*b22)  !!! CHECK !!!
      enddo

c ---- surface at j=jmx

      j=jmx
      do i=2,ii-1
       dphidxi(i) = 0.5*(phi(i+1,j)-phi(i-1,j))  !lagged derivative
      enddo
 
      do i=2,ii-1
       b12 = (zx(i,j)*ex(i,j) + zy(i,j)*ey(i,j))
       b22 = (ex(i,j)**2 + ey(i,j)**2)
       phi(i,j) = (-b12*dphidxi(i) + (2.0*phi(i,j-1)-0.5*phi(i,j-2))*b22)/(1.5*b22)  !!! CHECK !!!
      enddo
      
c **************************************************************
c ----- scrapcode for initializing velocity potential
c **************************************************************

      do j=1,jmx
      do i=1,imx
       phi(i,j) = Uinf*x(i,j)
      enddo
      enddo
%}
%% Mach number fix scrap code
%{
      subroutine props
      include "common.inc"

C Subroutine to calculate fluid properties

      xmcut2 = 2.0**2     ! user specified max Mach number
      xmtar2 = xmcut2*(1. + 0.5*(gamma-1.0)*xmach**2)   ! xmach = freestream Mach number
     c               /(1. + 0.5*(gamma-1.0)*xmcut2)     ! gamma = specific heat ratio
        
      do j= 1,jmax
      do i= 1,imax

       .
       .
       .

      u(i,j) = (zx(i,j)*dphidz + ex(i,j)*dphide)/xj(i,j)  ! initial calc. of U and V
      v(i,j) = (zy(i,j)*dphidz + ey(i,j)*dphide)/xj(i,j)

      vmagmax = sqrt(xmtar2*ainf**2)                      ! ainf = freestream sound speed
      vmag = sqrt(u(i,j)**2 + v(i,j)**2)

      u(i,j) = min(vmagmax,vmag)/(vmag+1e-6)*u(i,j)       ! correction to limit Mach number
      v(i,j) = min(vmagmax,vmag)/(vmag+1e-6)*v(i,j)

       .
       .
       .

      enddo
      enddo

%}
