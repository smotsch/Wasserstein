% This code computes the Wasserstein distance between two density functions
% Author: Shane Lubold, shane.lubold@asu.edu

% clear all; close all; clc
%% Insert the two functions to compute the Wasserstein distance between
% them.

%---------------------------------------------------
% Parameter Wasserstein distance
p = 2;                                  % Power of the Wasserstein distance

xPoints = [2 2 3 3];
yPoints = [1 4];
[rhoMacro,rhoMicro] = StochasticODEModel;
%xPoints = rhoMicro; yPoints = rhoMacro;


  xPoints = sort(rhoMicro);
  yPoints = sort(rhoMacro);
%xPoints = sort(xPoints);
%yPoints = sort(yPoints);


nX = length(xPoints);
nY = length(yPoints);

%---------------------------------------------------

% Define the cumulative functions, F, G
% if x_i, y_i set of points
 F = [0:nX]/nX;
 x = [xPoints xPoints(end)];
 G = [0:nY]/nY;
 y = [xPoints xPoints(end)];


% xDisc, rhoMicro
% F = cumsum(rhoMicro)*dx;
% G = cumsum(rhoMacro)*dx;
% x = [a-dx intX b+dx];
% y = [a-dx intX b+dx];


% F = [0 F 1];
% G = [0 G];

pCombined = sort(union(F, G)); 

% Construct F^-1, G^-1
FINV = interp1(F,x,pCombined,'left');
GINV = interp1(G,y,pCombined,'left');

%% Wasserstein distance using the trapezoidal method 
%WD = nthroot( trapz(pCombined,abs(FINV-GINV).^p) , p);
diff_FG = FINV(1:(end-1)) - GINV(1:(end-1));
WD_p = sum( abs(diff_FG).^p.*diff(pCombined) );
WD = nthroot( WD_p , p);
disp(['The Wasserstein distance is ',num2str(WD),'.']);
figure; plot(1:numel(rhoMacro),F,1:numel(rhoMicro),F); grid on;
title('\rho Macro vs. \rho Micro'); legend('Macro','Micro');

