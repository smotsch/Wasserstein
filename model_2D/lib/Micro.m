     function rho = Micro(x,y,D,dt,T,N,mu1,mu2,var1,var2)
% 
% Solve the KPP equation as a particle system:
%    dX_t = ?dB_t
%  + 'birth process' with intensity 1-?.
% 

%  function [rho] = KPP_Micro(rhoIC,x,D,dt,T,N,seedNum)
% 
% Solve the KPP equation as a particle system:
%    dX_t = ?dB_t
%  + 'birth process' with intensity 1-?.

%     x= linspace(-2,2,10);
%     y = x;
%     D= 100;
%     T = 1000;
%     dt = .1;
%     N = 10^3;
%     mu1 = 1;
%     mu2 = 0;
%     var1 = 1;
%     var2 = 2;

    [xmesh,ymesh] = meshgrid(x,y);
    rho = zeros(numel(x),numel(y),floor(T/dt + .5));
    
    X = [mu1 + sqrt(var1)*randn(N,1), mu2 + sqrt(var2)*randn(N,1)];
    rho(:,:,1) = hist3(X,[numel(x),numel(y)]);
   
    for j = 2:floor(T/dt + .5)
%         tic
        %  Move
        M = length(X);
        X = X + sqrt(2*D)*sqrt(dt)*randn(M,2); 
            
        % Histogram
        rho(:,:,j) = hist3(X,[numel(x),numel(y)]); 
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
    
%     figure;
%     for j = 1:floor(T/dt + .5)
%         surf(rho(:,:,j)); grid on;
%         xlabel('x'); ylabel('y'); zlabel('f(x,y)');
%         pause(.1);
%     end
%    

 end

