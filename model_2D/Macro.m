%function rho = Macro(x,y,D,dt,T)

% This code solves the equation: du/dt = D(d^2u/dx^2 + d^2u/dy^2) + u(1-u) 
% using a Crank-Nicolson scheme. Take dx = dy for simplicity. This code
% calls blktridiag, (http://www.mathworks.com/matlabcentral/fileexchange/
% 10603--block--tri-diagonal-matrices), which quickly constructs a
% tridiagional matrix from submatrices. The resulting equation is
% Au^{n+1} = Bu^{n} + f(u^{n}) so u^{n+1} = inv(A)*(Bu^{n}+dt*f(u))

% Author: Shane Lubold (shane.lubold@asu.edu).

clc; clear all; close all;


    % Model
    f = @(y) .1*y.*(1-y);

    x = linspace(-5,5,50);
    y = linspace(-6,6,60);
    D =2;
    dt = .1;
    T = 10;
    intT = 0:dt:T;
    [X1,X2] = meshgrid(x,y);
    X1_tp = X1';
    X2_tp = X2';
    rhoIC = mvnpdf([X1_tp(:) X2_tp(:)],[0, 0],[2,1]);

    
    % init
    dx = x(2)-x(1);
    ld = D*dt/dx^2;
    nX = length(x);
    nY = length(y);
    nT = floor(T/dt + .5); 
    tau = D*(dt/2*dx^2);
    
    rhoICmat =  reshape(rhoIC,nX,nY,1);
    
    % Initial Condition
    %rho(:,1) = reshape(rhoIC,nX*nY,1);
    rho = rhoIC;
    
    % Crank-Nicolson Method
    diagMatrixA =(1+4*tau)*diag(ones(1,nX)) - tau*diag(ones(1,nX-1),1) - tau*diag(ones(1,nX-1),-1);
    SubMatrixA = -tau*diag(ones(1,nX));
    A = full(blktridiag(diagMatrixA,SubMatrixA,SubMatrixA,nY)); %Construct "A"
    Ainv = inv(A); % Construct "A^{-1}"
    
    diagMatrixB =(1-4*tau)*diag(ones(1,nX)) + tau*diag(ones(1,nX-1),1) + tau*diag(ones(1,nX-1),-1);
    SubMatrixB = tau*diag(ones(1,nX));
    B = full(blktridiag(diagMatrixB,SubMatrixB,SubMatrixB,nY)); %Construct "B"

    rhoFINAL(:,:,1) = reshape(rhoIC,nX,nY,1);
    for j = 2:nT
        rho(:,j) = Ainv*(B*rho(:,j-1) + dt*f(rho(:,j-1)));
        rhoFINAL(:,:,j) = reshape(rho(:,j),nX,nY,1);
        surf(X1,X2,rhoFINAL(:,:,j)'); title('Solution');
        grid on; zlim([0 1.5]); xlabel('x'); 
        ylabel('y'); zlabel('f(x,y)'); 
        axis('equal')
        %a = gca; a.Box = 'on'; a.BoxStyle = 'full';
        colorbar;
        pause(.1);
    end

    %end

