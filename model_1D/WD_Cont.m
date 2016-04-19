function [WD,G,xG] = WD_Cont(f,g,x)
% 
% Compute the Wasserstein distance between two density distributions f and g defined on the same interval.
% 

    % A) compute the cdf
    [F,xF] = cdf_clean(x,f);
    [G,xG] = cdf_clean(x,g);
    
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

 
 