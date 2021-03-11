%% CFD Nathan Holmeswork 4

% Load Data
% Have 81 i and j vals
% 81 x 81 meshgrid
% imx and jmx not included in data
data = load('grid-poisson-2015');
%plot(data(2:end,1),data(2:end,2))
x = reshape(data(2:end,1),[81,81]);
y = reshape(data(2:end,2),[81,81]);
%plot(x,y)
imx = data(1,1);
jmx = data(1,2);

% Problem 2: phi
phi = 10.*x - 5.*y;
%% Create Computational Space Matricies
%{
for j = 1:jmx
    for i = 1:imx
        zx(i,j) = 0.5*(y(i,j+1) - y(i,j-1));
        zy(i,j) = -0.5*(x(i,j+1) - x(i,j-1));
        ex(i,j) = -0.5*(y(i+1,j) - y(i-1,j));
        ey(i,j) = 0.5*(x(i+1,j) - x(i-1,j));
        xj(i,j) = zx(i,j)*ey(i,j) - ex(i,j)*zy(i,j);
    end
end
%}

for j = 1:jmx
    for i = 1:imx
        % Eta terms
        if(i == 1)
            ex(i,j) = -(-1.5*y(i,j) + 2*y(i+1,j) - 0.5* y(i+2,j));
            ey(i,j) = (-1.5*x(i,j) + 2*x(i+1,j) - 0.5*x(i+2,j));
            
            % Problem 2
            phiz(i,j) = -1.5*phi(i,j) + 2*phi(i+1,j) - 0.5* phi(i+2,j);
        elseif (i == imx)
            ex(i,j) = -(1.5*y(i,j) - 2*y(i-1,j) + 0.5* y(i-2,j));
            ey(i,j) = (1.5*x(i,j) - 2*x(i-1,j) + 0.5*x(i-2,j));
            
            % Problem 2
            phiz(i,j) = 1.5*phi(i,j) - 2*phi(i-1,j) + 0.5*phi(i-2,j);
        else
            ex(i,j) = -0.5*(y(i+1,j) - y(i-1,j));
            ey(i,j) = 0.5*(x(i+1,j) - x(i-1,j));
            
            %Problem 2
            phiz(i,j) = 0.5*(phi(i+1,j) - phi(i-1,j));
        end
        
        % Zeta terms
        if(j == 1)
            zx(i,j) = -(-1.5*y(i,j) + 2*y(i,j+1) - 0.5* y(i,j+2));
            zy(i,j) = (-1.5*x(i,j) + 2*x(i,j+1) - 0.5*x(i,j+2));
            
            % Problem 2
            phie(i,j) = -(-1.5*phi(i,j) + 2*phi(i,j+1) - 0.5*phi(i,j+2));
        elseif (j == imx)
            zx(i,j) = -(1.5*y(i,j) - 2*y(i,j-1) + 0.5* y(i,j-2));
            zy(i,j) = (1.5*x(i,j) - 2*x(i,j-1) + 0.5*x(i,j-2));
            
            % Problem 2
            phie(i,j) = -(1.5*phi(i,j) - 2*phi(i,j-1) + 0.5*phi(i,j-2));
        else
            zx(i,j) = -0.5*(y(i,j+1) - y(i,j-1));
            zy(i,j) = 0.5*(x(i,j+1) - x(i,j-1));
            
            % Problem 2
            phie(i,j) = -0.5*(phi(i,j+1) - phi(i,j-1));
        end
        
        %xj = zx.*ey - ex.*zy;
    end
    % Calculation of xj array
    %xj(i,j) = zx(i,j)*ey(i,j) - ex(i,j)*zy(i,j);
    %u = (zx(i,j).*phiz(i,j) + ex(i,j).*phie(i,j))./xj(i,j);
    %v = (zy(i,j).*phiz(i,j) + ey(i,j).*phie(i,j))./xj(i,j);
end
xj = zx.*ey - ex.*zy; % calc xj array

% Problem 2
% u and v calculation: divide by 1/J (which is xj)
u = (zx.*phiz + ex.*phie)./xj;
v = (zy.*phiz + ey.*phie)./xj;
%% Plotting
%{
figure(1)
colormap('hsv')
hold on
contourf(x,y,ey,100,'edgecolor','none') % 100 contours at each val of ex
hold off
%}

% Problem 1
figure(1) % problem one contourf plots
colormap('jet')

%subplot(r,c,position number)
subplot(3,2,1) % ex
contourf(x,y,ex,100,'edgecolor','none')
title('\eta_x')
colorbar

subplot(3,2,2) % ey
contourf(x,y,ey,100,'edgecolor','none')
title('\eta_y')
colorbar

subplot(3,2,3) % zx
contourf(x,y,zx,100,'edgecolor','none')
title('\zeta_x')
colorbar

subplot(3,2,4) % zy
contourf(x,y,zy,100,'edgecolor','none')
title('\zeta_y')
colorbar

subplot(3,2,5) % xj
contourf(x,y,xj,100,'edgecolor','none')
title('xj')
colorbar

% Problem
figure(2) % u and v velocity components
colormap('jet');
subplot(1,2,1) % u component, subtract initial value (10)
contourf(x,y, u-10 ,100,'edgecolor','none')
title('u velocity component')
colorbar
subplot(1,2,2) % v component, subtract initial value (-5)
contourf(x,y, v-(-5) ,100,'edgecolor','none')
title('v velocity component')
colorbar

%% Problem 2
% Initilize phi grid like x and y using equation given
% Put FDA bit of phi in the loops
% Calculate u and v using equation format given

