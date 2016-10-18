function [densityMic,intR,int_rhoR] = radial_density(x,y,rhoMicro,rhoMacro,Part)

% Parameters
dr = .05;
intR = 0:dr:5;

% Initialize mean vectors
x_barMac = zeros(1,2); 
x_barMic = zeros(1,2);
rhoTempX = zeros(length(x),length(y));
rhoTempY = zeros(length(x),length(y));

for k = 1:length(rhoMicro(1,:,1))
rhoTempX(:,k) = x.*rhoMacro(k,:);
end

for k = 1:length(rhoMacro(:,1,1))    
rhoTempY(k,:) = y.*rhoMacro(:,k)';
end

% Compute x_bar
x_barMac(1) = 1/(trapz(x,trapz(y,rhoMacro))) * trapz(x,trapz(y,rhoTempX));
x_barMac(2) = 1/(trapz(x,trapz(y,rhoMacro))) * trapz(x,trapz(y,rhoTempY));

f = @(x,y) (x.^2+y.^2<1);
Mesh_x_barMac1 = x_barMac(1)*ones(length(x),length(x));
Mesh_x_barMac2 = x_barMac(2) *ones(length(y),length(y));
indicator_DiscR = @(x,y,R) ((x - Mesh_x_barMac1).^2 + (y - Mesh_x_barMac2).^2 < R^2);

% init
[X,Y] = meshgrid(x,y);
%A = f(X,Y);
A = exp(-X.^2-Y.^2);
nR = length(intR);
int_rhoR = zeros(1,nR);


% loop
for k=1:nR
    % init
    r = intR(k);
    B = indicator_DiscR(X,Y,r);
    % the value
    int_rhoR(k) = trapz(x,trapz(y,A.*B));
end

% plot
% intR_mid = (dr/2):dr:(2-dr/2);
% rhoR = diff(int_rhoR)./(2*pi*intR_mid*dr);
% grid on;

% plot(intR_mid,rhoR, intR_mid,exp(-intR_mid.^2))
% xlabel('r')

% Micro
x_barMic(1) = mean(Part(:,1));
x_barMic(2) = mean(Part(:,2));

for k=1:nR
    % init
    r = intR(k);
    densityMic(k) = sum((Part(:,1)-x_barMic(1)).^2 + (Part(:,2) - x_barMic(2)).^2  < r^2);
end
end




