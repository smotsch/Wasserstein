 function WD = WD_discreet_cont(X_points,intX,rhoMacro)
% 
% Compute the Wasserstein distance between a set of points {X}_i and a density distribution g.
% 
    % A) compute the cdf
    % A.1) cdf of the distribution of points
    nX = length(X_points);
    F = [0:(nX-1)]/(nX-1);              
    xF = sort(X_points);
    % A.2) as usual for G

    [xG,G] = cdf_clean(intX,rhoMacro);
    
    % B) Quantile
    z = sort(union(F, G));
    iMin = sum(z<1e-14);
    if (iMin>0)
        z(1:iMin) = [];
    end
    iMax = sum(z>1-1e-14);
    if (iMax>0)
        z((end-iMax):end) = [];
    end
    % C) inverse F^-1, G^-1
    FINV = interp1(F,xF,z,'left');
    GINV = interp1(G,xG,z,'left');
        
    % D) Wasserstein distance using the trapezoidal method 
    p = 2;
    WD = nthroot( trapz(z,abs(FINV-GINV).^p) , p);
 end
