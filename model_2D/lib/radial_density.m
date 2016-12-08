function FinalANS = radial_density(x,y,rhoMicro,rhoMacro,Part,intR)



% Initialize mean vectors
x_barMac = zeros(1,2); 
x_barMic = zeros(1,2);
rhoTempX = zeros(length(x),length(y));
rhoTempY = zeros(length(x),length(y));

for k = 1:length(rhoMicro(:,1))
rhoTempX(:,k) = x.*rhoMacro(k,:);
end

for k = 1:length(rhoMacro(1,:))    
rhoTempY(k,:) = y.*rhoMacro(:,k)';
end

% Compute x_bar
x_barMac(1) = 1/(trapz(x,trapz(y,rhoMacro))) * trapz(x,trapz(y,rhoTempX));
x_barMac(2) = 1/(trapz(x,trapz(y,rhoMacro))) * trapz(x,trapz(y,rhoTempY));

Mesh_x_barMac1 = x_barMac(1)*ones(length(x),length(x));
Mesh_x_barMac2 = x_barMac(2) *ones(length(y),length(y));
indicator_DiscR = @(x,y,R) ((x - Mesh_x_barMac1).^2 + (y - Mesh_x_barMac2).^2 < R^2);

% init
[X,Y] = meshgrid(x,y);
nR = length(intR);

% Macro
for k=1:nR
    % init
    r = intR(k);
    B = indicator_DiscR(X,Y,r);
    % the value
    densityMac(k) = trapz(x,trapz(y,rhoMacro.*B));
    
end

% Micro
x_barMic(1) = mean(Part(:,1));
x_barMic(2) = mean(Part(:,2));

for k=1:nR
    % init
    r = intR(k);
    densityMic(k) = sum((Part(:,1)-x_barMic(1)).^2 + (Part(:,2) - x_barMic(2)).^2  < r^2);
        
end

% Normalize Micro "CDF"
densityMic = densityMic./densityMic(end);



% plot(intR,densityMic/max(densityMic))
% title('Distribution of Particles')
% grid on;
% 
% diffDensity = diff(densityMic);
% plot(intR(2:end), diffDensity);
% title('Difference in Distribution of Particles')
% grid on;

%


FinalANS = [densityMic; densityMac]';
end
