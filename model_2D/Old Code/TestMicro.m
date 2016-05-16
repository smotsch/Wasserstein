 function rho = TwoDMicro(x,y, xInit, yInit,D,dt,T,N)
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

    
    [xInit, yInit] = AdjustSize(xInit', yInit');   
    X = interp2(rand(N,2),[xInit ; yInit]);  
    rho(:,:,1) =  hist3(X,[length(x),length(y)]);
    rho(:,:,1) = rho(:,:,1)./(sum(sum(rho(:,:,1))));

    X1 = X(:,1);
    X2 = X(:,2);
    
    for j = 2:floor(T/dt + .5)
        %  Move
        
        X1 = X1 + sqrt(2*D)*sqrt(dt)*randn(length(X1),1); 
        X2 = X2 + sqrt(2*D)*sqrt(dt)*randn(length(X2),1);     
        [X1, X2] = AdjustSize(X1,X2);
            
        % Histogram
        rho(:,:,j) = hist3([X1 X2],[length(x) length(y)]);
        rho(:,:,j) = rho(:,:,j)./(sum(sum(rho(:,:,j))));   
        
        rho_at_X = interp2([X1 X2],rho(:,:,j)); % interpolate at x_i, coming from Brownian paths.
        
        %% 2) Birth/Death Process  
        X1 = BirthDeath(X1,rho_at_X(:,1),dt);
        X2 = BirthDeath(X2,rho_at_X(:,2),dt);
       
    end
 end