% function TwoDStochastic
% This function computes rhoMicro, rhoMacro for a given initial condition
% in two dimensions. It calls Micro.m and Macro.m. 

% Author: Shane Lubold (shane.lubold@asu.edu).
clc; clear all; close all;

    dx = .2;
    dy = .2;
    dt = .1;
    T = 1;
    t = 0:dt:T;
    x = -1:dx:1;
    y = -1:dy:1;
%     x= linspace(-5,5,7);
%     y = x;
    D = 1;
    N = 10^2;
    c = .005; % threshold for computing contours.
    
    % Initial Conditiion for rhoMacro
    mu1 = -1;
    mu2 = 3;
    var1 = 1;
    var2 = 3;
    
    % Call Micro, Macro
    rhoMicro = Micro(x,y,D,dt,T,N); 
    rhoMacro = Macro(x,y,D,dt,T,mu1,mu2,var1,var2);
    
    % Compute Velocities
    for j = 1:floor(T/dt + .5)
        distanceMacro(j) = ComputeVelocity(x,y,rhoMacro(:,:,j),c);
    end
    
    Fit = polyfit(1:size(rhoMacro,3),distanceMacro./dt,1);
    figure; plot(Fit);
    xlabel('Time'); ylabel('Velocity');
    
    %Compute the WD
    [WD,Diff1,Diff2] = TwoDOptimizationCode(x,y,rhoMicro(:,:,end),...
    rhoMacro(:,:,end))
    
  

%     subplot(2,2,1); surf(x,y,rhoMicro(:,:,1)'); colorbar;
%     xlabel('x'); ylabel('y'); title('Micro, Initial');
%     subplot(2,2,2); surf(x,y,rhoMacro(:,:,1)'); colorbar;
%     xlabel('x'); ylabel('y'); title('Macro, Initial')
%     subplot(2,2,3); surf(x,y,rhoMicro(:,:,end)');colorbar;
%     xlabel('x'); ylabel('y'); title('Micro, End')
%     subplot(2,2,4); surf(x,y,rhoMacro(:,:,end)');colorbar;
%     xlabel('x'); ylabel('y'); title('Macro, End')
    
% end