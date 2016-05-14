% function rho = TwoDMicro(rhoIC,x,D,dt,T,N,seedNum)
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

    % init
    x = -4:.2:4;
    y = -4:.2:4;
    T = 1;
    dt = .1;
    D = 1;
    dx = x(2)-x(1);
    dy = y(2) - y(1);
    nX = length(x);
    nY = length(y);
    nT = floor(T/dt + .5); 
    rho = zeros(nX,nY,nT+1);
    sigma = sqrt(2*D);
    N = nX;
    
      xInit = normpdf(x,0,2);
      yInit = normpdf(y,0,1);
      X = interp2(rand(N,2),[xInit ; yInit]);
      
      Histogram =  hist3(X,[nX,nY]);
      rho(:,:,1) = (Histogram)./(sum(sum(Histogram)));

    %-------------------------------------------%
    %---            Big loop                 ---%
    %-------------------------------------------%
    tic
    
    for k = 1:nT
        % A) Move
        X1 = X(:,1);
        X2 = X(:,2);
        
        if k ~= 1
            X1 = cell2mat(X(:,1));
            X2 = cell2mat(X(:,2));
        end
        
        M = length(X1);
        K = length(X2);
        
        X1 = X1 + sigma*sqrt(dt)*randn(M,1);
        X2 = X2 + sigma*sqrt(dt)*randn(K,1);
        
        % Make them the size for purposes of hist3, interp2
        if k~=1  
            lengthDesired = max(length(X1),length(X2));
            if length(X1) == lengthDesired
                increase = length(X1) - length(X2);
                adjust = NaN(1,increase);
                X2 = [X2; adjust'];
            end
            if length(X2) == lengthDesired
                increase = length(X2) - length(X1);
                adjust = NaN(1,increase);
                X1 = [X1; adjust'];
            end
        end
        
        % B) Birth/death
        rho(:,:,k+1) = hist3([X1 X2],[nX, nY],'DisplayStyle','Tiles');
        rho(:,:,k+1) = rho(:,:,k+1)./(sum(sum(rho(:,:,k+1)))); % Normalize
        
        
        rho_at_X = interp2([X1 X2],rho(:,:,k+1)); % interpolate at x_i, coming from Brownian paths.
        %% ∂_t ρ = ρ(1-ρ)
        %% 1) density
        % Then make them back to their original size.        
        
        oneMrho1 = (1-rho_at_X(:,1));               % 1-M at X(k)       
        oneMrho2 = (1-rho_at_X(:,2));
        %% 2) rand      
        M = length(X1);
        K = length(X2);
        
        coin_M1 = rand(M,1);
        coin_M2 = rand(K,1);
        
        indexKill1     = logical( (oneMrho1<0).*(coin_M1<dt*abs(oneMrho1)) );
        indexKill2     = logical( (oneMrho2<0).*(coin_M2<dt*abs(oneMrho2)) );
        
        indexDivision1 = logical( (oneMrho1>0).*(coin_M1<dt*abs(oneMrho1)) );
        indexDivision2 = logical( (oneMrho2>0).*(coin_M2<dt*abs(oneMrho2)) );
        
        %% 3) update
        if (~isempty(indexDivision1))  
            % First do column 1
           X_new1 = X1(indexDivision1);
        end
        if (~isempty(indexDivision1))
            % First do column 1
            X_new2 = X2(indexDivision2);
        end
        
        if (~isempty(indexKill1))
            X1(indexKill1) = [];
        end                 
        if (~isempty(indexKill2))
            X2(indexKill2) = [];
        end
        
        if (~isempty(indexDivision1))
            X1 = [X1; X_new1];
        end
        if (~isempty(indexDivision2))
            X2 = [X2; X_new2];
        end   
    X = {X1 X2};
    end
    
  figure;
  for k = 1:nT
    surf(rho(:,:,k)); grid on; colorbar;
    pause(.3);
  end