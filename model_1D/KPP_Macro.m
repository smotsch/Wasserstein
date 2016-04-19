function rho = KPP_Macro(rhoIC,x,D,dt,T,methodNum)
% 
% Solve the KPP equation:
%    ∂_t ρ = DΔ_x ρ + ρ(1-ρ)
% We use Dirichlet boundary condition (ρ = 0)
% 
% Input:
%    . rhoIC, x  : initial condition ρ_0(x) on the interval x
%    . D         : coefficient of diffusion
%    . dt, T     : time
%    . methodNum : 1-> explicit, 2-> Cranck-Nicholson
% Output:
%    . rho       : the numerical solution ρ(x,t)
%    

    % Model
    f = @(y) y.*(1-y);
    
    % init
    dx = x(2)-x(1);
    ld = D*dt/dx^2;
    nT = floor(T(end)/dt + .5); 
    nX = length(x);
    rho = zeros(nX,nT+1);

    % Initial Condition
    rho(:,1) = rhoIC;
    
    % Init method
    if (methodNum == 1)
        % numerical trick: explicit method
        intM = 2:(nX-1);		% points on the Middle
        intR = 3:nX;		        % points on the Right
        intL = 1:(nX-2);		% points on the Left
    else
        % Cranck-Nicholson
        tao = D*(dt/2)/dx^2;

        % Create A & B
        A = -tao.*diag(ones(1,nX-1),-1) +...
            (1+2.*tao).*diag(ones(1,nX))...
            -tao.*diag(ones(1,nX-1),1);
        
        B = tao*diag(ones(1,nX-1),-1) + ...
            (1-2*tao).*diag(ones(1,nX)) + ... 
            tao*diag(ones(1,nX-1),1);
        % Neuman BC
        A(end,1) = -tao; A(1,end) = -tao;
        B(end,1) =  tao; B(1,end) =  tao;
        % ρ_n+1 = A^(-1)*B ρ_n
        Ainv_B = inv(A)*B;
    end    
    
    
    %-------------------------------------------%
    %---            Big loop                 ---%
    %-------------------------------------------%

    if (methodNum == 1)
        % explict method
        for k = 1:nT
            rho(:,k+1) = (1-2*ld)*rho(:,k) + ld*([rho(2:end,k); rho(1,k)] + [rho(end,k); rho(1:(end-1),k)]) ...
                + dt*f(rho(:,k));
        end    
    else
        % Cranck-Nicholson
        for k = 1:nT
            % Strang splitting
            % ----------------
            % A) ∂_t u = f(u) over Δt/2
            rho1 = rho(:,k) + dt/2*f(rho(:,k));
            % B) ∂_t u = ∂_xx u over Δt
            rho2 = Ainv_B*rho1;
            % C) ∂_t u = f(u) over Δt/2
            rho(:,k+1) = rho2 + dt/2*f(rho2);
        end
    end
    %-------------------------------------------%
    %-------------------------------------------%
    
end



% remark:
%--------
%  with BC ρ=0 at the boundary
%  ---------------------------
%  intM = 2:(nX-1);		% points on the Middle
%  intR = 3:nX;			% points on the Right
%  intL = 1:(nX-2);		% points on the Left
%   
%  for k = 1:nT
%      rho(intM,k+1) = (1-2*ld)*rho(intM,k) + ld*(rho(intL,k) + rho(intR,k)) ...
%          + dt*f(rho(intM,k));
%  end    
%  
