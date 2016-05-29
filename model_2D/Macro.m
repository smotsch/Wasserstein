function rhoFINAL = Macro(x,y,D,dt,T,mu1,mu2,var1,var2,shouldPlot)

% This code solves the equation: du/dt = D(d^2u/dx^2 + d^2u/dy^2) + u(1-u) 
% using a Crank-Nicolson scheme. Take dx = dy for simplicity. This code
% calls blktridiag, (http://www.mathworks.com/matlabcentral/fileexchange/
% 10603--block--tri-diagonal-matrices), which quickly constructs a
% tridiagional matrix from submatrices. The resulting equation is
% Au^{n+1} = Bu^{n} + f(u^{n}) so u^{n+1} = inv(A)*(Bu^{n}+dt*f(u))

% Author: Shane Lubold (shane.lubold@asu.edu).

    % Model
    f = @(y) y.*(1-y);

    [X1,X2] = meshgrid(x,y);
    X1_tp = X1';
    X2_tp = X2';
    rhoIC = mvnpdf([X1_tp(:) X2_tp(:)],[mu1, mu2],diag([var1, var2]));
        
    % init
    dx = x(2)-x(1);
    nX = length(x);
    nY = length(y);
    nT = floor(T/dt + .5); 
    tau = D*dt/(2*dx^2);
    
    % Initial Condition
    rho = rhoIC;
    rhoFINAL = zeros(nX,nY,nT+1);
    rhoFINAL(:,:,1) = reshape(rhoIC,nX,nY,1);
    
    % Crank-Nicolson Method
    addpath('lib')
    diagMatrixA =(1+4*tau)*diag(ones(1,nX)) - tau*diag(ones(1,nX-1),1) - tau*diag(ones(1,nX-1),-1);
    SubMatrixA = -tau*diag(ones(1,nX));
    A = full( blktridiag(diagMatrixA,SubMatrixA,SubMatrixA,nY) ); 
    Ainv = inv(A); 
    diagMatrixB =(1-4*tau)*diag(ones(1,nX)) + tau*diag(ones(1,nX-1),1) + tau*diag(ones(1,nX-1),-1);
    SubMatrixB = tau*diag(ones(1,nX));
    B = full( blktridiag(diagMatrixB,SubMatrixB,SubMatrixB,nY) ); 


    % loop time
    for j = 2:nT
        rho = Ainv*(B*rho + dt*f(rho));
        rhoFINAL(:,:,j) = reshape(rho,nX,nY,1);
        % plot
        if (shouldPlot)
            surf(X1,X2,rhoFINAL(:,:,j)'); 
            title(['Solution at t=',num2str(j*dt,'%1.2f')]);
            grid on; xlabel('x'); 
            zlim([0 1.5]);
            ylabel('y'); zlabel('f(x,y)'); 
            axis('equal')
            %a = gca; a.Box = 'on'; a.BoxStyle = 'full';
            colorbar;
            pause(.1);
        end
    end

end

