  function [xPos1] = CalculateVelocity(f,intX,xLoc1,xLoc2)
% This code computes the velocity of a function.
% fix a a value of Y as the "reference point". Author: Shane Lubold
% shane.lubold@asu.edu

yStar = .1;                            %threshold
xM = find(intX==0);
fTemp = f(xM:end);           %cut and use only positive x values
[~, loc1] = min(abs(fTemp-yStar)); 
intXTemp = intX(xM:end);
xPos1 = intXTemp(loc1);
end
  