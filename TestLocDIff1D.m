% This code solves the diffusion equation in one dimension using a
% diffusion coefficient that depends on the the physical location. 

x = linspace(-3,3,10);
t = linspace(0,1,.1);
D = linspace(0,5,numel(x));
rho = zeros(numel(x), numel(t));
deltaX = x(2)-x(1);

rho(:,1) = normpdf(x,0,1);


for i = 2:numel(t)
    for j = 2:numel(x)
       rho(
    end
end