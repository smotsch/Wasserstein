% function TwoDStochastic
% This function computes rhoMicro, rhoMacro for a given initial condition
% in two dimensions. It calls Micro.m and Macro.m. 

% Author: Shane Lubold (shane.lubold@asu.edu).
clc; clear all; close all;

    dx = .2;
    dy = .2;
    dt = .01;
    T = 1;
    t = 0:dt:T;
    x = -4:dx:3;
    y = -4:dy:4;
    D = 50;
    N = 1000;
    c = .005; % threshold for computing contours.
    
    % Initial Conditiion for rhoMacro
    mu1 = 0;
    mu2 = 0;
    var1 = 2;
    var2 = 2;
    
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
  
    
 



%     subplot(2,2,1); surf(x,y,rhoMicro(:,:,1)'); colorbar;
%     xlabel('x'); ylabel('y'); title('Micro, Initial');
%     subplot(2,2,2); surf(x,y,rhoMacro(:,:,1)'); colorbar;
%     xlabel('x'); ylabel('y'); title('Macro, Initial')
%     subplot(2,2,3); surf(x,y,rhoMicro(:,:,end)');colorbar;
%     xlabel('x'); ylabel('y'); title('Micro, End')
%     subplot(2,2,4); surf(x,y,rhoMacro(:,:,end)');colorbar;
%     xlabel('x'); ylabel('y'); title('Macro, End')
    
% end