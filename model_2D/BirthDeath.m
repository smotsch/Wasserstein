function X = BirthDeath(X,Interpolation,dt)

% This function takes in vectors X, Interpolation and a scalar dt. 
% It creates new cells and kills cells according to some probability.

% Author: Shane Lubold (shane.lubold@asu.edu)

oneMrho = 1-Interpolation;               % 1-M at X(k) 
coin_M = randn(length(X),2);

indexKill = logical( (oneMrho<0).*(coin_M<dt*abs(oneMrho)) );
indexDivision = logical( (oneMrho>0).*(coin_M<dt*abs(oneMrho)) );

X_new = X(indexDivision);
X(indexKill) = [];

X = [X; X_new];
X(isnan(X)) = [];
end

