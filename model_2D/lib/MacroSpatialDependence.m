
% This code solves the equation: du/dt = D(d^2u/dx^2 + d^2u/dy^2) + alpha*u(1-u)
% where D and alpha are spatially dependent

% Author: Shane Lubold (shane.lubold@asu.edu).
    clc; clear all; close all;
  
    % Take in Blood Vessel Data
    RawData = imread('RawData.jpg');
    RawDataScaled = RawData(:,:,1)/max(max(RawData(:,:,1)));
    [numelX, numelY] = size(RawDataScaled);
  
    T = 1;
    dt = .1;
    x = linspace(-1,1,numelX);
    y = linspace(-1,1,numelY);
    mu1 = 0;
    mu2 = 0;
    var1 = .2;
    var2 = .2;    
    dx = x(2)-x(1);
    dy = y(2)-y(1);
    nX = length(x);
    nY = length(y);
    nT = floor(T/dt + .5);
   
    [X1,X2] = meshgrid(x,y);
    X1_tp = X1';
    X2_tp = X2';
    rho = zeros(length(x),length(y),nT);

    rhoIC = reshape(mvnpdf([X1_tp(:) X2_tp(:)],[mu1, mu2],[var1, var2]),nX,nY);

    % If we want Dirichlet BC (boundaries @ 0), use the following code
    BC = 0; 
    rho(:,:,1) = rhoIC;
    rho(1,:,:) = BC; rho(:,1,:) = BC; rho(end,:,:) = BC; rho(:,end,:) = BC; %hold the boundaries at zero.
   
   % Build the D and alpha
   a = 1; 
   b = -3;
   c = 1; 
   d = 1;
   
   alpha = @(z) a+b*z;
   D = @(z) c+d*z;
 
    tic
    for k = 1:nT %k for time
        for i = 2:length(x)-1 %i for x
            for j = 2:length(y)-1 % j for y
                rho(i,j,k+1) = rho(i,j,k) + D(RawDataScaled(i+1, j))*(dt/dx^2)*(rho(i+1,j,k) - rho(i,j,k)) ...
                    + D(RawDataScaled(i-1, j))*(dt/dx^2)*(rho(i+1,j,k) - rho(i-1,j,k)) ... 
                + D(RawDataScaled(i, j+1))*(dt/dy^2)*(rho(i,j+1,k) - rho(i,j,k)) +... 
                D(RawDataScaled(i, j-1))*(dt/dy^2)*(rho(i,j+1,k) - rho(i,j-1,k)) + dt * alpha(RawDataScaled(i, j)) * ...
                rho(i,j,k)*(1-rho(i,j,k));     
            end
        end
    end
    toc
    
    
    
