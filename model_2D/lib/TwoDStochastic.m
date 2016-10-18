% function TwoDStochastic
% This function computes rhoMicro, rhoMacro for a given initial condition
% in two dimensions. It calls Micro.m and Macro.m. 

% Author: Shane Lubold (shane.lubold@asu.edu).
clc; clear all; close all;

%% Parameters
D = 1;                                % diffusion coefficients
%% Parameters micro
N = 10^2;                               % number of cells (initially)
Mass_rhoIC = 2;
seedNum = 1;
% domain
dx = .1;
dy = .1;
%x = -5:dx:5;
x = linspace(-5,5,2^3);
y = x;
%y = -5:dy:5;
% time
dt = .01;
T =   1;
t = 0:dt:T;

% Initial Conditiion
mu1 = 0;
mu2 = 0;
var1 = 1;
var2 = 1;
shouldPlot = 0;
r = linspace(1.5,40,100);
VelMac = zeros(length(t),length(r));
VelMic = zeros(length(t),length(r));




[rhoMicro, Part] = Micro(x,y,D,dt,T,N,mu1,mu2,var1,var2);
[rhoDir,rho] = Macro(x,y,D,dt,T,mu1,mu2,var1,var2,true);
for l = 1:length(t)
   VelMac(l,:) = ComputeVelocity(x,y,rho(:,:,l),r); 
end
 

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
                                  
[A, b, C] = OptimMatrices(rhoMicro(:,:,q),rho(:,:,q),x,y);
 % Compare to solution from linprog
 tic
 [xval, fval(q),exitflag] = linprog(C,[],[],A,b,zeros(size(C)),[]);
 LinProgCode = toc;
   
 end
 xval = reshape(xval,length(x)^2, length(y)^2);
 

 
% end