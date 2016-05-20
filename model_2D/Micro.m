%  function rho = Micro(x,y, xInit, yInit,D,dt,T,N)
% 
% Solve the KPP equation as a particle system:
%    dX_t = σdB_t
%  + 'birth process' with intensity 1-ρ.
% 

% function [rho] = KPP_Micro(rhoIC,x,D,dt,T,N,seedNum)
% 
% Solve the KPP equation as a particle system:
%    dX_t = σdB_t
%  + 'birth process' with intensity 1-ρ.
% 

    clc; clear all; close all;
    dx = .1;
    dy = .1;
    dt = .1;
    x = -3:dx:3;
    y = -3:dy:3;
    N = 10;
    T = 0:dt:50;
    nT = floor(T(end)/dt + .5);
    dt = .1;
    D = 1;
    [xmesh,ymesh] = meshgrid(x,y);
    rho = zeros(numel(x),numel(y),nT);
    
    X = [rand(N,1), rand(N,1)];
    rho(:,:,1) = hist3(X,[numel(x),numel(y)]);
   
    for j = 2:nT
        %  Move
        M = length(X);
        X = X + sqrt(2*D)*sqrt(dt)*randn(M,2); 
            
        % Histogram
        rho(:,:,j) = hist3(X,[numel(x),numel(y)]);   
        rho_at_X = interp2(xmesh,ymesh,rho(:,:,j),X(:,1),X(:,2)); % interpolate at x_i, coming from Brownian paths.
        
        %% 2) Birth/Death Process  
        oneMrho = 1-rho_at_X;               % 1-M at X(k) 
        coin_M = rand(M,1);
        
        indexKill = logical( (oneMrho<0).*(coin_M<dt*abs(oneMrho))); 
        indexDivision = logical( (oneMrho>0).*(coin_M<dt*abs(oneMrho)) );

        if (length(indexDivision)>0) % New Cells
            X_new = X(indexDivision,:);
        end
         
        if (length(indexKill)>0)
            X(indexKill,:) = [];
        end

        if (length(indexDivision)>0)
            X = [X; X_new];
        end
    
    end
    
    figure;
    for j = 1:nT
        surf(rho(:,:,j)); grid on;
        xlabel('x'); ylabel('y'); zlabel('f(x,y)');
        pause(.1);
    end
   
%  end