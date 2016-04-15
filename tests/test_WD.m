%% test Gaussian
%%--------------

addpath('../model_1D')

% init
dx = .001;
x = -20:dx:20;

sig1 = 1; m1 = 3;
sig2 = 2; m2 = 2;

f = exp(-(x-m1).^2/(2*sig1^2))/sqrt(2*pi*sig1^2);
g = exp(-(x-m2).^2/(2*sig2^2))/sqrt(2*pi*sig2^2);


% Numerical estimation
WD_Cont(f,g,x)
% theoretical value for Gaussian 
sqrt( (m1-m2)^2+(sig1-sig2)^2 )

% Numerical estimation
X = m1 + sig1*randn(20000,1);
WD_discreet_cont(X,g,x)