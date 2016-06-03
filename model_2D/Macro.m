<<<<<<< HEAD
% function rhoFINAL = Macro(x,y,D,dt,T,mu1,mu2,var1,var2)
=======
function rhoFINAL = Macro(x,y,D,dt,T,mu1,mu2,var1,var2,shouldPlot)
>>>>>>> origin/master

% This code solves the equation: du/dt = D(d^2u/dx^2 + d^2u/dy^2) + u(1-u)
% with two approaches. The first uses finite differencing two solve under
% both Dirichlet and Neumann BC. The second uses a Crank-Nicolson under 
% an unknown BC. This code
% calls blktridiag, (http://www.mathworks.com/matlabcentral/fileexchange/
% 10603--block--tri-diagonal-matrices), which quickly constructs a
% tridiagional matrix from submatrices. The resulting equation is
% Au^{n+1} = Bu^{n} + f(u^{n}) so u^{n+1} = inv(A)*(Bu^{n}+dt*f(u))

% Author: Shane Lubold (shane.lubold@asu.edu).
    clc; clear all; close all;
    f = @(u) u.*(1-u);    
    dx = 1.5;
    dy = 1.5;
    x = -5:dx:5;
    y = -5:dy:5;
    mu1 = 0;
    mu2 = 0;
    var1 = .2;
    var2 = .2;
    nX = length(x);
    nY = length(y);
    T = 20;
    dt = 10^-2;
    D = 1;
    nT = floor(T/dt + .5);
    [X1,X2] = meshgrid(x,y);
    X1_tp = X1';
    X2_tp = X2';
<<<<<<< HEAD
    rhoIC = reshape(mvnpdf([X1_tp(:) X2_tp(:)],[mu1, mu2],[var1, var2]),nX,nY);

    % If we want Dirichlet BC (boundaries @ 0), use the following code
    BC = 0; 
    rhoDir(:,:,1) = rhoIC;
    rhoDir(1,:,:) = BC; rhoDir(:,1,:) = BC; rhoDir(end,:,:) = BC; rhoDir(:,end,:) = BC; %hold the boundaries at zero. 
 
    tic
    for k = 1:nT-1 %k for time
        for i = 2:length(x)-1 %i for x
            for j = 2:length(y)-1 % j for y
                rhoDir(i,j,k+1) = rhoDir(i,j,k) + (dt/dx^2)*(rhoDir(i+1,j,k) -2*rhoDir(i,j,k) + rhoDir(i-1,j,k)) ...
                + (dt/dy^2)*(rhoDir(i,j+1,k) -2*rhoDir(i,j,k) + rhoDir(i,j-1,k)) + dt*rhoDir(i,j,k)*(1-rhoDir(i,j,k));     
            end
        end
    end
    toc
    
%     surf(rho(:,:,end-1)); grid on %plot the resulting surface
%     rhoMacro = reshape(rho(:,:,end-1),numel(rho(:,:,end-1)),1); %reshape the vector to pass to other function. 

% If we want Neumann BC (derivative on boundary = 0), use the following code

tic
rho(:,:,1) = rhoIC;
for k = 1:nT-1
=======
    rhoIC = mvnpdf([X1_tp(:) X2_tp(:)],[mu1, mu2],diag([var1, var2]));
        
    % init
    dx = x(2)-x(1);
    nX = length(x);
    nY = length(y);
    nT = floor(T/dt + .5); 
    tau = D*dt/(2*dx^2);
>>>>>>> origin/master
    
   rho(1,1,k+1) =  rho(1,1,k)+ (dt/dx^2)*(rho(2,1,k)-rho(1,1,k)) + (dt/dy^2)*(rho(1,2,k)-rho(1,1,k));
   rho(end,end,k+1) =  rho(end,end,k)+ (dt/dx^2)*(-rho(end,end,k)+rho(end-1,end,k)) + (dt/dy^2)*(-rho(end,end,k)+rho(end,end-1,k));
   rho(end,1,k+1) =  rho(end,1,k)+ (dt/dx^2)*(-rho(end,1,k)+rho(end-1,1,k)) + (dt/dy^2)*(rho(end,2,k)-rho(end,1,k));
   rho(1,end,k+1) =  rho(1,end,k)+ (dt/dx^2)*(rho(2,end,k)-rho(1,end,k)) + (dt/dy^2)*(-rho(1,end,k)+rho(1,end-1,k));
   
   for j = 2:nX-1
       rho(1,j,k+1) = rho(1,j,k) + (dt/dx^2)*(rho(2,j,k)-rho(1,j,k)) +(dt/dy^2)*(rho(i,j+1,k)-2*rho(i,j,k)+rho(i,j-1,k)); 
       rho(end,j,k+1) = rho(end,j,k) + (dt/dx^2)*(-rho(end,j,k)+rho(end-1,j,k)) + (dt/dy^2)*(rho(i,j+1,k)-2*rho(i,j,k)+rho(i,j-1,k)); 
   end

   for i = 2:nY-1
       rho(i,1,k+1) = rho(i,1,k) + (dt/dx^2)*(rho(i+1,j,k)-2*rho(i,j,k)+rho(i-1,j,k)) + (dt/dy^2)*(rho(i,2,k)-rho(i,1,k));
       rho(i,end,k+1) = rho(i,end,k) + (dt/dx^2)*(rho(i+1,j,k)-2*rho(i,j,k)+rho(i-1,j,k)) + (dt/dy^2)*(-rho(i,end,k)+rho(i,end-1,k)); 
   end

   for i = 2:length(x)-1 %i for x
        for j = 2:length(y)-1 % j for y
                rho(i,j,k+1) = rho(i,j,k) + (dt/dx^2)*(rho(i+1,j,k) -2*rho(i,j,k) + rho(i-1,j,k)) ...
                + (dt/dy^2)*(rho(i,j+1,k) -2*rho(i,j,k) + rho(i,j-1,k)) + dt*rho(i,j,k)*(1-rho(i,j,k));     
        end
   end
end

toc



% %     
% 
%     
%     
    % Initial Condition
<<<<<<< HEAD
    rho2(:,1) = reshape(rhoIC, numel(rhoIC),1);
    alpha = (dt/(2*dx^2));
    beta = (dt/(2*dy^2));
    NumEl = numel(x);
    
    % Crank-Nicolson Method
    diagMatrixA =(1+2*alpha+2*beta)*diag(ones(1,NumEl)) - alpha*diag(ones(1,NumEl-1),1) - alpha*diag(ones(1,NumEl-1),-1);
    SubMatrixA = -beta*diag(ones(1,NumEl));
    A = full(blktridiag(diagMatrixA,SubMatrixA,SubMatrixA,NumEl)); 
    Ainv = inv(A); 

    diagMatrixB =(1-2*alpha-2*beta)*diag(ones(1,NumEl)) + alpha*diag(ones(1,NumEl-1),1) + alpha*diag(ones(1,NumEl-1),-1);
    SubMatrixB = beta*diag(ones(1,NumEl));
    B = full(blktridiag(diagMatrixB,SubMatrixB,SubMatrixB,NumEl)); 
    
    rhoFINAL(:,:,1) = reshape(rhoIC,NumEl,NumEl,1);
    tic
    for j = 2:nT
        rho2(:,j) = Ainv*(B*rho2(:,j-1) + dt*f(rho2(:,j-1)));
        rhoFINAL(:,:,j) = reshape(rho2(:,j),NumEl,NumEl,1);
               
%         surf(X1,X2,rhoFINAL(:,:,j)'); title('Solution');
%         grid on; xlabel('x'); 
%         zlim([0 1.5]);
%         ylabel('y'); zlabel('f(x,y)'); 
%         axis('equal')
%         a = gca; a.Box = 'on'; a.BoxStyle = 'full';
%         colorbar;
%         pause(.1);
=======
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
>>>>>>> origin/master
    end
    toc

% end



    
