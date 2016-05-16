function TwoDStochastic
% This function computes rhoMicro, rhoMacro for a given initial condition
% in two dimensions. It calls Micro.m and Macro.m. 

% Author: Shane Lubold (shane.lubold@asu.edu).
clc; clear all; close all;

    dx = .2;
    dy = .2;
    dt = .1;
    T = 1;
    x = -4:dx:3;
    y = -4:dy:4;
    D = 1;
    N = 10000;
    
    
    mu1 = 1;
    mu2 = 0;
    var1 = 1;
    var2 = 2;
    
    xInit = normpdf(x,mu1,var1);
    yInit = normpdf(y,mu2,var2);
    
    % Call Micro, Macro
    rhoMicro = Micro(x,y,xInit,yInit,D,dt,T,N); 
    rhoMacro = Macro(x,y,D,dt,T,mu1,mu2,var1,var2);
    
    subplot(2,2,1); surf(x,y,rhoMicro(:,:,1)'); colorbar;
    xlabel('x'); ylabel('y'); title('Micro, Initial');
    subplot(2,2,2); surf(x,y,rhoMacro(:,:,1)'); colorbar;
    xlabel('x'); ylabel('y'); title('Macro, Initial')
    subplot(2,2,3); surf(x,y,rhoMicro(:,:,end)');colorbar;
    xlabel('x'); ylabel('y'); title('Micro, End')
    subplot(2,2,4); surf(x,y,rhoMacro(:,:,end)');colorbar;
    xlabel('x'); ylabel('y'); title('Macro, End')
    
end