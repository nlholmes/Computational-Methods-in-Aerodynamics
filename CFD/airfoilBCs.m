function [phi,x,y] = airfoilBCs(phi,ex,ey,zx,zy,imax,jmax,x,y,farfield)
%   Deals with tangency boundary condition for CFD project
    %farfield = phi.*x; % square
    %phi(:,jmax) = farfield(:,jmax);
    % Tangency boundary condition
    j=1; % only do j=1,not j=jmax
    dphidxi = zeros(imax-1,1);
    for i = 2:imax-1
        dphidxi(i) = 0.5*(phi(i+1,j)-phi(i-1,j));  %!lagged derivative
    end

    for i = 2:imax-1 % 1 to 2, no change
        %dphidx(i) = 0.5*(phi(i+1,j)-phi(i-1,j));
        b12 = (zx(i,j)*ex(i,j) + zy(i,j)*ey(i,j));
        b22 = (ex(i,j)^2 + ey(i,j)^2);
        phi(i,j) = (-b12*dphidxi(i) + (2.0*phi(i,j+1)-0.5*phi(i,j+2))*b22)/(1.5*b22);  %!!! CHECK !!!
    end
    
    % Recent BC notes
    % Tangency, Far field, wake cut
    % Wake cut, doesnt seem to do anything here
    
    
    x(1,:) = x(imax-2,:);
    y(1,:) = y(imax-2,:); % -2 to -1
    x(imax,:) = x(3,:);
    y(imax,:) = y(3,:);
    
    
    % Far field --> perhaps do first like in channel (caused it to work)
    phi(:,jmax) = farfield(:,jmax);
    %phi(imax,:) = farfield(1,:);
    
end

