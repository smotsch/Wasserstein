function distance = ComputeVelocity(x,y,rhoMacro,c)

[Y,I] = min(abs(rhoMacro(:)-c));
[row, col] = ind2sub(size(rhoMacro),I);
distance = sqrt(x(row)^2 + y(col)^2);

end