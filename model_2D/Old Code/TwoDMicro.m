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
        
        % Make them the size for purposes of hist3
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
        rho(:,:,k+1) = (hist3([X1 X2],[nX, nY]));
        rho(:,:,k+1) = rho(:,:,k+1)./(sum(sum(rho(:,:,k+1)))); % Normalize
        
        rho_at_X = interp2([X1 X2],rho(:,:,k+1)); % interpolate at x_i, coming from Brownian paths.
        %% ∂_t ρ = ρ(1-ρ)
        %% 1) density
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
toc
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



