function [F,xF]=cdf_clean(x,f)
% Compute a strictly increasing CDF:
%  F(x) = ∫_-∞^x f(y) dy
% The set [x] is a uniform grid.

% init
dx = x(2) - x(1);
% A) Basic cdf
%-------------
F = cumsum(f)*dx;
xF = x;
% B) Fine tuning
%---------------
% B.1) From  0 à 1
xF = [x(1)-dx xF];
F = F/F(end);

F = [0; F];
% B.2) Cut value above 1 (due to numerical errors)
F = min(F,1);
% B.3) Remove multiple valu e
ind_multiple_F = find(diff(F)<1e-20);
F(ind_multiple_F) = [];
xF(ind_multiple_F) = [];
end

