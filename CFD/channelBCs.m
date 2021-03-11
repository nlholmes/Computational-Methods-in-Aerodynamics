function [phi] = channelBCs(phi,ex,ey,zx,zy,imax,jmax,farfield)
%   Deals with tangency boundary condition for CFD project 2


    % Recent BC notes
    % Tangency, Far field, wake cut
    % Far field
    phi(imax,:) = farfield(imax,:); % right side
    phi(1,:) = farfield(1,:); % left
    phi(:,jmax) = farfield(:,jmax); % top
    phi(1,:) = farfield(1,:); % bottom
    %{
    c **************************************************************
    c ----- scrapcode for applying surface tangency bcs
    c ----- set V-dot-n = 0 =>b12*dphidxi + b22*dphideta = 0.0
    c **************************************************************
    %}
    %c ---- surface at j=1

    % Tangency boundary condition
    j=1;
    for i = 2:imax-1
        dphidx(i) = 0.5*(phi(i+1,j)-phi(i-1,j));  %!lagged derivative
    end

    for i = 2:imax-1
        %dphidx(i) = 0.5*(phi(i+1,j)-phi(i-1,j));
        b12 = (zx(i,j)*ex(i,j) + zy(i,j)*ey(i,j));
        b22 = (ex(i,j)^2 + ey(i,j)^2);
        phi(i,j) = (-b12*dphidx(i) - (2.0*phi(i,j+1)-0.5*phi(i,j+2))*b22)/(-1.5*b22);  %!!! CHECK !!!
    end

    %c ---- surface at j=jmx

    j=jmax;
    for i = 2:imax-1
        dphidx(i) = 0.5*(phi(i+1,j)-phi(i-1,j));  %!lagged derivative
    end

    for i = 2:imax-1
        b12 = (zx(i,j)*ex(i,j) + zy(i,j)*ey(i,j));
        b22 = (ex(i,j)^2 + ey(i,j)^2);
        % removed neg before b12
        % added neg sign before 2*phi (from plus)
        % -b12
        % removed neg sign before phi and denom
        phi(i,j) = (-b12*dphidx(i) + (2.0*phi(i,j-1)-0.5*phi(i,j-2))*b22)/(1.5*b22);  %!!! CHECK !!!
    end
    
    
    
    % first row x first row y farfield initial phi
    %{      
    c **************************************************************
    c ----- scrapcode for initializing velocity potential
    c **************************************************************

          do j=1,jmx
          do i=1,imx
           phi(i,j) = Uinf*x(i,j)
          enddo
          enddo
    %}
end

