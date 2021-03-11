%% CFD Project 1
% LWE

%% Preallocation/Input section

%%%%%% Run this section before running any other sections %%%%%%

% Inputs
length = 1; % from 0.25 to 0.85
dist = 0.6; % distance to travel
deltaX = 0.01;
wave_speed = 1; % variable c
nu = 1; % cfl/vN number, 0.4, 1.0, 1.3
time = dist/wave_speed; % total time to run

% Dependent inputs
deltaT = nu*deltaX/wave_speed;
nsteps = time/deltaT;

% Number of points
imx = length/deltaX + 1;

% Preallocation
x = ones(1,imx);
u = zeros(1,imx);
utmp = zeros(1,imx);
% Generate 1D mesh
for i = 1:imx
    x(i) = deltaX*(i-1);
end

% Initialize solution
for i = 1:imx
    % This if sequence is slightly redundant
    if x(i) < 0.1
        u(i) = 0;
    elseif x(i) >= 0.1 && x(i) <= 0.25
        u(i) = 1;
    else
        u(i) = 0;
    end
end
u = u'; % transpose for implicit method

%% Explicit backward spacial difference, first order in time
% Loop over time
for n = 1:10
    for i = 2:imx-1
       utmp(i) = u(i) - nu*(u(i) - u(i-1));
    end
    u(:,:) = utmp(:,:);
    plot(x(:,:),u(:,:));
    hold on
    pause(0.1)
    ylim([0,5])
end

%% Explicit central spacial difference, first order in time
% From class --> for any nu, will never be stable (von neumann analysis)
% Loop over time
for n = 1:5
    for i = 2:imx-1
       utmp(i) = u(i) - nu/2*(u(i+1) - u(i-1));
       %utmp(i) = u(i+1) - nu*(u(i+1) - u(i-1)); % this works but thats weird, nu = 1
        %utmp(i) = u(i-1); % this is why
    end
    u(:,:) = utmp(:,:);
    %u(1,2:imx-1) = utmp(1,2:imx-1);
    plot(x(:,:),u(:,:));
    
    hold on
    pause(0.1)
end

%% Implicit central spacial difference, first order in time
% Set Boundary Values
BV1 = 0; % Start boundary value for (u_i)^(n+1)
BV2 = 0; % End boundary value

% Preallocate
a = zeros(1,imx-2);
d = a;
c = a;
b = a;
X = a;
% Create tridiagonal matrix vectors
%{
for i = 1:imx-2
    a(i) = -nu/2;
    d(i) = 1;
    c(i) = -a(i);
    for i = 2:imx-1
        b(i-1) = u(i);
    end
end
% Assign boundary values to b vector
b(1) = b(1) + nu/2*BV1;
b(imx-2) = b(imx-2) + nu/2*BV2;
%}
for n = 1:nsteps
    for i = 1:imx-2
        a(i) = -nu/2;
        d(i) = 1;
        c(i) = -a(i);
        for i = 2:imx-1
            b(i-1) = u(i);
        end
    end
    % Assign boundary values to b vector
    %b(1) = b(1) + nu/2*BV1;
    %b(imx-2) = b(imx-2) + nu/2*BV2;
    b(1) = 0;
    b(imx-2) = 0;
    
    % Solve tridiagonal matrix using thomas algorithm
    for i = 2:(imx-2)
        w = a(i)/d(i-1);
        d(i) = d(i) - w*c(i-1);
        b(i) = b(i) - w*b(i-1);
    end
    X(imx-2) = b(imx-2)/d(imx-2);
    for i = (imx-3):-1:1
        %X(i+1) = b(i+1)/d(i+1);
        X(i) = (b(i) - c(i)*X(i+1))/d(i);
    end
    
    for a = 2:imx-3
    u(a) = X(a);
    end
    
    u(1) = 0;
    u(imx-2) = 0;
    plot(x,u)
    xlim([0,0.6]);
    hold on
    pause(.1)
end
%% Explicit Lax-Wendroff, second order in time and space
% Loop over time
for n = 1:nsteps
    for i = 2:imx-1
       utmp(i) = u(i) - 1/2*nu*(u(i+1) - u(i-1)) + 1/2*nu^2*(u(i+1) - 2*u(i) + u(i-1));
    end
    u(:,:) = utmp(:,:);
    plot(x(:,:),u(:,:));
    hold on
    pause(0.1)
end

%% Implicit method using backslash solve
tridiag = full(gallery('tridiag',imx,-nu/2,1,nu/2));
for i = 1:nsteps
    utmp = tridiag\u;
    u = utmp;
    plot(x,u)
    hold on
    pause(0.1)
end