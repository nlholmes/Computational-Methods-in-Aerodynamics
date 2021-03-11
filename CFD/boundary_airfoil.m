function [phi, x, y] = boundary_airfoil(x, y, phi, ex, ey, zx, zy, imax, jmax, farfield)
%applies conditions at the boundary 
%  

j = 1; % define phi for j = 1 
    for i = 2:imax-1
        dphidxi(i) = 0.5*(phi(i+1,j)-phi(i-1,j)); % lagged derivative
    end
    for i = 2:imax-1
        b12 = zx(i,j)*ex(i,j) + zy(i,j)*ey(i,j);
        b22 = (ex(i,j).^2 + ey(i,j).^2);
        phi(i,j) = (-b12*dphidxi(i) - (2.0*phi(i,j+1)-0.5*phi(i,j+2))*b22)/(-1.5*b22);  % CHECK !!!
    end

    % Wake Cut condition
    x(1,:) = x(imax-2, :); 
    y(1,:) = y(imax-2, :); % potential error, not -2, changed and improved bounds of contour
    x(imax, :) = x(3,:); 
    y(imax,:) = y(3,:); 
    
    % Farfield condition
    phi(:,jmax) = farfield(:,jmax);
end