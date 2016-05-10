% function rho = Macro(x,y,D,dt,T)

% This code solves the equation: du/dt = D(d^2u/dx^2 + d^2u/dy^2) + u(1-u) 
% using a Crank-Nicolson scheme. Take dx = dy.


    % Model
    f = @(y) y.*(1-y);
    
    % init
    dx = x(2)-x(1);
    ld = D*dt/dx^2;
    nX = length(x);
    nY = length(y);
    nT = floor(T/dt + .5); 
    rho = zeros(nX,nY,nT+1);
    tau = D*(dt/2*dx^2);
    % Initial Condition
%     rho(:,1) = rhoIC;
    
    % Crank-Nicolson Method
    A =(1+2*tau)*diag(ones(1,nX)) - tau*diag(ones(1,nX-1),1) - tau*diag(ones(1,nX-1),-1);
    B = -tau*diag(ones(1,nX));
    C = [B A B];
    FINAL(1:nX,1:nX) = A;
    FINAL(nX+1:2*nX,nX+1:2*nX) = B;
    for j =2:nX^2-2
        for i = 1:nX^2
            if
        FINAL()    
        end
    end
    
   
    
    
    
    
% end
