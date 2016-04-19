% Authors: Shane Lubold (shane.lubold@asu.edu) & Sebastien Motsch
% (smotsch@asu.edu).

% This code takes two Gaussians and computes the Wasserstein distance
% between them using the following methodology.
% Solve the optimal transport problem between two sets of points in R^2:
%    {X_i}_i   ?  {Y_j}_j
% The two sets can have different size (nX and nY), each X particle weights 1/nX 
% and the Y particle 1/nY (the 'total weight' is the same).
%  
% The goal is to find the optimal 'transport plan' (x_ij) where:
%     x_ij: part of particle X_i move to X_j
% It minimizes the cost (L^2):
%     ?_ij x_ij |Y_j-X_i|^2
% and satisfies the constraints
%     ?_j x_ij = 1/nX
%     ?_i x_ij = 1/nY
% or in matrix notation:
%    A X = b
% with A=[1 1 1 ...; ...] and b=[1/nX .. 1/nX 1/nY .. 1/nY].

% The algorithm used the simplex method:
%  . Phase 1: find a 'basic' solution
%      -> Change AX=b into A~X=b~ with A~ containing columns of the form (0..0 1 0 ..0).
%      -> We deduce a basic solution X = (0 ..0 * 0 .. 0 * .. ..0)
%  . Phase 2: solve the minimization problem with the constraints A~X=b~

% A) initial configuration
%-------------------------

    clear all; 
    clc; close all;

    mu_f    = [0 0];
    Sigma_f = [1 .3; .3 1];
    mu_g    = [2 2];
    Sigma_g = [.25 .3; .3 1];
    
    [X1,X2] = meshgrid(linspace(-5,5,7)',linspace(-5,5,7)');
    XMESH = [X1(:) X2(:)];
    f = mvnpdf(XMESH,mu_f,Sigma_f);
    g = mvnpdf(XMESH,mu_g,Sigma_g);
    
    FSUM = sum(sum(f)); GSUM = sum(sum(g));
    f = f/FSUM; g = g/GSUM;
    
    figure; 
    surf(X1,X2,reshape(f,7,7)); hold on
    surf(X1,X2,reshape(g,7,7)); hold off;
    
    Xpos = XMESH;
    Ypos = Xpos;
    nX = length(Xpos);
    nY = length(Ypos);
    
    % B) constraint:
    %---------------
    % We are lookin for a transition matrix  '(x_ij)' with
    %   x_ij piece of i going to j
    % '(x_ij)' is a nX×nY matrix.  
    % Thus,
    % ?_i x_ij = 1/nX   -> each Xpos particle 'weight' 1/nX
    % ?_j x_ij = 1/nY   -> each Ypos particle 'weight' 1/nY
    % We transform the matrix '(a_ij)' into a vector 'X'
    
    b = [f; g];
    A = zeros(nX+nY,nX*nY);
    for k=1:nX
        A(k, (k-1)*nY+[1:nY]) = 1;
    end
    for k=1:nY
        for k2=1:nX
            A(k+nX, k + (k2-1)*nY) = 1;
        end
    end
    % We have an extra constraint (last equality is automatically satisfied)
     A(end,:) = [];
     b(end)   = [];

    % C) cost matrix
    %---------------
    % |Y_k'-X_k|^2
    nX = length(Xpos(:,1));
    nY = length(Ypos(:,1));
    matrix_C = zeros(nX,nY);
    for i=1:nX
        for j=1:nY
            matrix_C(i,j) = sqrt( (Ypos(j,1)-Xpos(i,1))^2 + (Ypos(j,2)-Xpos(i,2))^2 );
        end
    end
    C = matrix_C(:);

    %---------------------------------------------------------------%
    %-----------------         Resolution            ---------------%
    %---------------------------------------------------------------%

    %----- Phase 1 -----%
    %-------------------%
    fprintf('---------------------------------------\n')
    fprintf('Phase 1: find a basic feasible solution\n')
    fprintf('---------------------------------------\n')
%     if (nX==nY)
%         % simple basic feasible solution: 
%         %    Xi is sent to Yi
%         %   and Xn sent to no-one else
%         %  -> better to use the Hungarian algorithm in this case
%         %     http://cstheory.stackexchange.com/questions/2707/justification-for-the-hungarian-method-kuhn-munkres
%         tic
%         matrix_X = diag(ones(nX,1));
%         X        = 1/nX*matrix_X(:);
%         A_phase2 = A;
%         b_phase2 = b;
%         for i=1:(nX-1)
%             A_phase2(i+nX,:) = A_phase2(i+nX,:) - A_phase2(i,:);
%             b_phase2(i+nX,:) = b_phase2(i+nX,:) - b_phase2(i,:);
%         end
%         for i=1:(nX-1)
%             A_phase2(nX,:) = A_phase2(nX,:) - A_phase2(nX+i,:);
%             b_phase2(nX,:) = b_phase2(nX,:) - b_phase2(nX+i,:);
%         end    
%         bv         = zeros(2*nX-1,1);
%         bv(1:nX)   = [1 1+(nX+1)*(1:(nX-1))];
%         bv((nX+1):end) = nX^2-[(nX-1):(-1):1];
%         timePhase1=toc
%     else
        [A_phase2,b_phase2,X,bv] = MinimizeTableau(A,b,-sum(A),-sum(b));
%     end
    % reduce gradient 'r'
    r = C';
    for i=1:length(bv)
        % find associated basic variable
        j = bv(i);
        if (j>0)
            % update cost
            %  ?_x_N r = c_N' - c_b'?(A_b^-1 A_n)
            %          = c'   - c_b'?(A_b^-1 A)
            r = r - C(j)*A_phase2(i,:);
        end
    end

    %----- Phase 2 -----%
    %-------------------%
    fprintf('---------------------------------------\n')
    fprintf('Phase 2: find the optimal solution\n')
    fprintf('---------------------------------------\n')
    initial_cost = C'*X;
    [A_opt,b_opt,X_opt,bv_opt] = MinimizeTableau(A_phase2,b_phase2,r,initial_cost,bv);
    final_cost = C'*X_opt

    %---------------------------------------------------------------%
    %---------------------------------------------------------------%

    % plot the result
    matrix_X_opt = reshape(X_opt,nY,nX)';
    clf
    hold on
    plot(Xpos(:,1),Xpos(:,2),'bo',Ypos(:,1),Ypos(:,2),'rx')
    for i=1:nX
        for j=1:nY
            if matrix_X_opt(i,j)>1d-10
                quiver([Xpos(i,1),0],[Xpos(i,2),0], ...
                       [Ypos(j,1)-Xpos(i,1),0],[Ypos(j,2)-Xpos(i,2),0],'autoScale','off')
            end
        end
    end
    hold off
    
    % Check to see if the constraints are being satisfied
    onesMatrixf = ones(length(Xpos),1);
    onesMatrixg = ones(length(Ypos),1);
    
    diff1 = f - matrix_X_opt*onesMatrixf;
    diff2 = g - (onesMatrixg'*matrix_X_opt)';
    
    % Calculate the true value of W_2
    diff_mean = sqrt((mu_g(2)-mu_f(2) + (mu_g(1)-mu_f(1) ) ) );
    True_WD = diff_mean + trace(Sigma_f) + trace(Sigma_g) - 2*trace((Sigma_g.^(1/2)*Sigma_g*Sigma_g.^(1/2))^(1/2))
    error = abs(True_WD - final_cost)
    
   