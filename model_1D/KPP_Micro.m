  function [rho,X] = KPP_Micro(rhoIC,x,D,dt,T,N,seedNum)
% 
% Solve the KPP equation as a particle system:
%    dX_t = σdB_t
%  + 'birth process' with intensity 1-ρ.
% 

    % init
    dx = x(2)-x(1);
    nX = length(x);
    nT = floor(T(end)/dt + .5); 
    rho = zeros(nX,nT+1);
    sigma = sqrt(2*D);

    % Initial Condition
    %   -> compute the cdf of ρ_0 to generate the points
    F = cumsum(rhoIC)*dx;
    Mass_rhoIC = F(end);
    F = F/F(end);
%     rand(seedNum);
%     rand(seedNum);
    X = interp1(F,x,rand(N,1));
    rho(:,1) = hist(X,x)*Mass_rhoIC/N/dx;

    %-------------------------------------------%
    %---            Big loop                 ---%
    %-------------------------------------------%
    
    for k = 1:nT
        % A) Move
        M = length(X);
        X = X + sigma*sqrt(dt)*randn(M,1);
        % B) Birth/death
        rho(:,k+1) = hist(X,x)*Mass_rhoIC/N/dx;
        rho_at_X = interp1(x,rho(:,k+1),X); % interpolate at x_i, coming from Brownian paths.
        %% ∂_t ρ = ρ(1-ρ)
        %% 1) density
        oneMrho = (1-rho_at_X);               % 1-M at X(k)
        %% 2) rand
        coin_M = rand(M,1);
        indexKill     = logical( (oneMrho<0).*(coin_M<dt*abs(oneMrho)) );
        indexDivision = logical( (oneMrho>0).*(coin_M<dt*abs(oneMrho)) );
        %% 3) update
        if (length(indexDivision)>0)
            X_new = X(indexDivision,:);
        end
        if (length(indexKill)>0)
            X(indexKill,:) = [];
        end
        if (length(indexDivision)>0)
            X = [X; X_new];
        end
    end
    %-------------------------------------------%
    %-------------------------------------------%
 end

