function f_int = ComputeVelocity(x,y,rho,r)

xRestrict = zeros(length(x),1); yRestrict = zeros(length(y),1);

for i = 1:length(r)
for k = 1:length(y)
    for j = 1:length(x)
        if (x(j).^2+y(k).^2 <= r(i).^2)
            xRestrict(j) = 1; yRestrict(k) = 1;
        end
    end    
end

f_res = rho(logical(xRestrict),logical(yRestrict),1);
f_int(i) = 1/(2*pi*r(i)) * trapz(x(find(xRestrict)),trapz(y(find(yRestrict)),f_res));

end


end