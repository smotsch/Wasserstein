function [rho,X] = Micro(x,y,D,dt,T,N,mu1,mu2,var1,var2)

% 
% Solve the KPP equation as a particle system:
%    dX_t = ?dB_t
%  + 'birth process' with intensity 1-?.
% 
    [X1,X2] = meshgrid(x,y);
    [xmesh,ymesh] = meshgrid(x,y);
    rho = zeros(numel(x),numel(y),floor(T/dt + .5));
    
    X = [mu1 + sqrt(var1)*randn(N,1), mu2 + sqrt(var2)*randn(N,1)];
    %rho(:,:,1) = reshape(mvnpdf([X1_tp(:) X2_tp(:)],[mu1, mu2],[var1, var2]),nX,nY); 
    rho(:,:,1) = hist3(X,[numel(x),numel(y)]);
   
    for j = 2:floor(T/dt + .5)
%         tic
        %  Move
        M = length(X);
        X = X + sqrt(2*D)*sqrt(dt)*randn(M,2); 
            
        % Histogram
        %rho(:,:,j) = hist3(X,[numel(x),numel(y)]); 
        [bandwidth, density, x_vec, y_vec] = kde2d(X,numel(x));
        rho(:,:,j) = density;
        [Xtemp1, Xtemp2] = meshgrid(X(:,1), X(:,2));
        rho_at_X = interp2(xmesh,ymesh,rho(:,:,j),Xtemp1(:,1),Xtemp2(:,1)); % interpolate at x_i, coming from Brownian paths.
        
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
%         toc
    
    end

 end

