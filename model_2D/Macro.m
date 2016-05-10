function rho = Macro(x,y,D,dt,T)

% This code solves the equation: du/dt = D(d^2u/dx^2 + d^2u/dy^2) + u(1-u) 
% using a Crank-Nicolson scheme. Take dx = dy for simplicity. This code
% calls blktridiag, (http://www.mathworks.com/matlabcentral/fileexchange/
% 10603--block--tri-diagonal-matrices), which quickly constructs a
% tridiagional matrix from submatrices. The resulting equation is
% Au^{n+1} = Bu^{n} + f(u^{n}) so u^{n+1} = inv(A)*(Bu^{n}+dt*f(u))

% Author: Shane Lubold (shane.lubold@asu.edu).

clc; clear all; close all;


    % Model
    f = @(y) y.*(1-y);

    x = linspace(-5,5,50);
    y = x;
    T = linspace(0,10,50);
    D =1;
    dt = .1;
    [X1,X2] = meshgrid(x,y);
    rhoIC = mvnpdf([X1(:) X2(:)],[0, 0],[1,1]);
    
    % init
    dx = x(2)-x(1);
    ld = D*dt/dx^2;
    nX = length(x);
    nY = length(y);
    nT = floor(T(end)/dt + .5); 
    tau = D*(dt/2*dx^2);
    % Initial Condition
    rho(:,1) = reshape(rhoIC,nX*nY,1);
    
    % Crank-Nicolson Method
    diagMatrixA =(1+4*tau)*diag(ones(1,nX)) - tau*diag(ones(1,nX-1),1) - tau*diag(ones(1,nX-1),-1);
    SubMatrixA = -tau*diag(ones(1,nX));
    A = full(blktridiag(diagMatrixA,SubMatrixA,SubMatrixA,nX)); %Construct "A"
    Ainv = inv(A); % Construct "A^{-1}"
    
    diagMatrixB =(1-4*tau)*diag(ones(1,nX)) + tau*diag(ones(1,nX-1),1) + tau*diag(ones(1,nX-1),-1);
    SubMatrixB = tau*diag(ones(1,nX));
    B = full(blktridiag(diagMatrixB,SubMatrixB,SubMatrixB,nX)); %Construct "B"

    rhoFINAL(:,:,1) = reshape(rhoIC,nX,nY,1);
    for j = 2:length(T)
        rho(:,j) = Ainv*(B*rho(:,j-1) + dt/2*f(rho(:,j-1)));
        rhoFINAL(:,:,j) = reshape(rho(:,j),nX,nY,1);
        surf(X1,X2,rhoFINAL(:,:,j)); title('Solution');
        grid on; zlim([0 1.5]); xlabel('x'); 
        ylabel('y'); zlabel('f(x,y)'); 
        a = gca; a.Box = 'on'; a.BoxStyle = 'full';
        colorbar;
        pause(.1);
    end

end

