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

    % init
    x = -4:.2:4;
    y = -4:.2:4;
    T = 10;
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
    
    for k = 1:nT
        % A) Move
        M = length(X(:,1));
        N = length(X(:,2));
        
        X(:,1) = X(:,1) + sigma*sqrt(dt)*randn(M,1);
        X(:,2) = X(:,2) + sigma*sqrt(dt)*randn(N,1);
        
        % B) Birth/death
        rho(:,:,k+1) = hist3(X,[nX, nY]);
        rho_at_X = interp2(X,rho(:,:,k+1)); % interpolate at x_i, coming from Brownian paths.
        %% ∂_t ρ = ρ(1-ρ)
        %% 1) density
        oneMrho = (1-rho_at_X);               % 1-M at X(k)       
        %% 2) rand
        coin_M = rand(M,2);
        indexKill     = logical( (oneMrho<0).*(coin_M<dt*abs(oneMrho)) );
        indexDivision = logical( (oneMrho>0).*(coin_M<dt*abs(oneMrho)) );
        
        %% 3) update
        if (length(indexDivision)>0)  
            % First do column 1
           X_new1 = X(indexDivision(:,1));
           X_new2 = X(indexDivision(:,2));
        end
        if (length(indexKill)>0)
            X(indexKill) = [];
        end                                
        if (length(indexDivision)>0)
            X1 = [X(:,1); X_new1];
            X2 = [X(:,2); X_new2];
            X = [X1 X2];
        end
    end
    %-------------------------------------------%
    %-------------------------------------------%
% end
























% Create Parameters
T_final = 15;
x1 = -5; 
x2 = 5;
y1 = -5; 
y2 = 5; 
N = 100;
intX = linspace(x1,x2,150);
intY = linspace(y1,y2,150); 
dx = intX(2) - intX(1);
dy = intY(2) - intY(1); 
t = linspace(0,T_final,100);
dt = t(2) - t(1);
sigma = 1;
[xFinal, yFinal] = meshgrid(intX,intY);
rhoIC = reshape(mvnpdf([xFinal(:), yFinal(:)],[0,0],[1,1]),length(intX),length(intY));
seedNum = 10;

% X = zeros(numIt,length(intX),length(intY)); %initialize matrix
% xBrown = zeros(length(t),length(intX));
% yBrown = zeros(length(t),length(intY)); 
% xBrown(1,:) = 0; yBrown(1,:) = 0;
% 
% randNumber1 = randn(numIt,length(intX),length(intX));
% randNumber2 = randn(numIt,length(intX),length(intY));
% %% Compute Brownian components x_i, y_i separately. 
% for i = 2:numel(t)
%%     X = X+sigma*sqrt(dt)*randn(N,2);
%     xBrown(1:N,i) = xBrown(1:numIt,i-1) + sigma*sqrt(dt)*randNumber1(1:numIt,i-1);  
%     yBrown(1:N,i) = yBrown(1:numIt,i-1) + sigma*sqrt(dt)*randNumber2(1:numIt,i-1);
% end
% 
% X = [xBrown; yBrown];
% figure; plot(xBrown,yBrown); grid on; 
% xlabel('x'); ylabel('y'); title('Brownian Paths');
% TwoDHist = hist3([xBrown(:,end),yBrown(:,end)]);
% surf(TwoDHist); grid on;

    
    
    % Initial Condition
    %   -> compute the cdf of ρ_0 to generate the points
    F = cumsum(rhoIC)*dx;
    Mass_rhoIC = F(end);
    F = F/F(end);
    rand('seed',seedNum);
    randn('seed',seedNum);
    X = interp1(F,[intX intY],rand(numIt,2),'nearest');
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
            X_new = X(indexDivision);
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

% end


%     % Initial Condition
%     %   -> compute the cdf of ρ_0 to generate the points
%     xInit = normpdf(x,0,2);
%     yInit = normpdf(y,0,1);
%     FInit = xInit'*yInit; %f(x,y) = f(x)*f(y), Independent Variables.
%     figure; 
%     surf(x',y',FInit);
%     CDF = histogram2(xInit,yInit,nX,'Normalization', 'cdf'); % Find the CDF
%     CDF = CDF.Values;
%     figure;
%     surf(x',y',CDF); xlabel('x'); ylabel('y'); zlabel('F(x,y)');
% 
%     [gridx,gridy] = meshgrid(rand(N,1),rand(N,1)); % find interpolation mesh
%     [xmesh, ymesh] = meshgrid(-4:.1:4,-4:.1:4);
%     X = interp2(CDF,gridx,gridy,'nearest');



