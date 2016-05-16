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
    
    x = linspace(-2,2,10);
    y = linspace(-2,2,10);
    xInit = normpdf(x,0,1);
    yInit = normpdf(y,0,1);
    N = 25;
    T = 1;
    dt = .1;
    D =1 ;
       
    X = interp2(rand(N,2),[xInit ; yInit]);  
    rho(:,:,1) =  hist3(X,[length(x),length(y)]);
    rho(:,:,1) = rho(:,:,1)./(sum(sum(rho(:,:,1))));

    
    for j = 2:floor(T/dt + .5)
        %  Move
        M = length(X);
        X = X + sqrt(2*D)*sqrt(dt)*randn(M,2); 
            
        % Histogram
        rho(:,:,j) = hist3(X,[length(x) length(y)]);   
        rho_at_X = interp2(X,rho(:,:,j)); % interpolate at x_i, coming from Brownian paths.
        
        %% 2) Birth/Death Process  
        oneMrho = 1-rho_at_X;               % 1-M at X(k) 
        coin_M = rand(length(X),1);

        indexKill1 = logical( (oneMrho(:,1)<0).*(coin_M<dt*abs(oneMrho(:,1))));
        indexKill2 = logical( (oneMrho(:,2)<0).*(coin_M<dt*abs(oneMrho(:,2))));
        indexKill = logical(indexKill1 .* indexKill2);
        indexDivision1 = logical( (oneMrho(:,1)>0).*(coin_M<dt*abs(oneMrho(:,1))) );
        indexDivision2 = logical( (oneMrho(:,2)>0).*(coin_M<dt*abs(oneMrho(:,2))) );
        indexDivision = logical(indexDivision1 .* indexDivision2);
        
        X_new = X([indexDivision, indexDivision]);
        if (~isempty(find(indexKill)))
            X(indexKill, indexKill) = [];
        end

        X = [X; X_new'];
       
    end
%  end