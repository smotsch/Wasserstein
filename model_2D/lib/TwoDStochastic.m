% function TwoDStochastic
% This function computes rhoMicro, rhoMacro for a given initial condition
% in two dimensions. It calls Micro.m and Macro.m. 

% Author: Shane Lubold (shane.lubold@asu.edu).
clc; clear all; close all;

%% Parameters
D = .1;                                % diffusion coefficients
%% Parameters micro
N = 10^3;                               % number of cells (initially)
Mass_rhoIC = 2;
seedNum = 1;
% domain
dx = .5;
dy = .5;
x = -2:dx:2;
y = -2:dy:2;
% time
dt = .1;
T =   20;
t = 0:dt:T;

% Initial Conditiion
mu1 = 0;
mu2 = 0;
var1 = 1;
var2 = 1;
shouldPlot = 0;


rhoMicro = Micro(x,y,D,dt,T,N,mu1,mu2,var1,var2);
[rhoDir,rho] = Macro(x,y,D,dt,T,mu1,mu2,var1,var2,true);
 

for q = 1:length(t)
% Call Micro, Macro
%------------------

rhoMicro(:,:,q) = rhoMicro(:,:,q)/sum(sum(rhoMicro(:,:,q)));
rho(:,:,q) = rho(:,:,q)/sum(sum(rho(:,:,q)));

% Compute speed of traveling waves
%---------------------------------
% addpath('lib')
% c = .2; % threshold for computing contours.
% for j = 1:floor(T/dt + .5)
%     distanceMacro(j) = ComputeVelocity(x,y,rho(:,:,j),c);
% end
% Fit = polyfit(t(2:end),distanceMacro,1);
% figure; 
% plot(t(2:end),distanceMacro,t(2:end),polyval(Fit,t(2:end)))
% xlabel('Time'); ylabel('Velocity');

% break
    
%Compute the WD
%tic 
%[WD(q),Diff1,Diff2,C,A,b] = TwoDOptimizationCode(x,y,rhoMicro(:,:,q),...
 %                                     rho(:,:,q));
%SimplexCode = toc;

f = rho(:,:,q);
g = rhoMicro(:,:,q);

    [X1, X2] = meshgrid(x,y);
    XMESH = [X1(:) X2(:)];
    Xpos = XMESH;
    Ypos = Xpos;
    nX = length(Xpos);
    nY = length(Ypos);

    b = [f(:); g(:)];
    
    A = zeros(nX+nY,nX*nY);
    for k=1:nX
        A(k, (k-1)*nY+[1:nY]) = 1;
    end
    for k=1:nY
        for k2=1:nX
            A(k+nX, k + (k2-1)*nY) = 1;
        end
    end
    
%      % We have an extra constraint (last equality is automatically satisfied)
      A(end,:) = [];
      b(end)   = [];
      
          nX = length(Xpos(:,1));
    nY = length(Ypos(:,1));
    matrix_C = zeros(nX,nY);
    for i=1:nX
        for j=1:nY
            matrix_C(i,j) = sqrt( (Ypos(j,1)-Xpos(i,1))^2 + (Ypos(j,2)-Xpos(i,2))^2 );
        end
    end
    C = matrix_C(:);
      

                                    
 % Compare to solution from linprog
 tic
 [xval, fval(q),exitflag] = linprog(C,[],[],A,b,zeros(size(C)),[]);
 LinProgCode = toc;
   
 end
 xval = reshape(xval,length(x)^2, length(y)^2);
 


 
% end