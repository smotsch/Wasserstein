function [final_cost, maxDiff1, maxDiff2] = TwoDptimizationCode(x,y,f,g)


    [X1, X2] = meshgrid(x,y);
    XMESH = [X1(:) X2(:)];
    Xpos = XMESH;
    Ypos = Xpos;
    nX = length(Xpos);
    nY = length(Ypos);

    b = [f(:); g(:)];
    
    A = zeros(nX+nY,nX*nY);
    for k=1:nX
        A(k, (k-1)*nY+[1:nY]) = 1;
    end
    for k=1:nY
        for k2=1:nX
            A(k+nX, k + (k2-1)*nY) = 1;
        end
    end
    
%      % We have an extra constraint (last equality is automatically satisfied)
      A(end,:) = [];
      b(end)   = [];
     
    nX = length(Xpos(:,1));
    nY = length(Ypos(:,1));
    matrix_C = zeros(nX,nY);
    for i=1:nX
        for j=1:nY
            matrix_C(i,j) = sqrt( (Ypos(j,1)-Xpos(i,1))^2 + (Ypos(j,2)-Xpos(i,2))^2 );
        end
    end
    C = matrix_C(:);
    
    fprintf('---------------------------------------\n')
    fprintf('Phase 1: find a basic feasible solution\n')
    fprintf('---------------------------------------\n')
    [A_phase2,b_phase2,X,bv] = MinimizeTableau(A,b,-sum(A),-sum(b));
    r = C';
    
    for i=1:length(bv)
        % find associated basic variable
        j = bv(i);
        if (j>0)
            % update cost
            r = r - C(j)*A_phase2(i,:);
        end
    end
    
    fprintf('---------------------------------------\n')
    fprintf('Phase 2: find the optimal solution\n')
    fprintf('---------------------------------------\n')
    initial_cost = C'*X;
    [A_opt,b_opt,X_opt,bv_opt] = MinimizeTableau(A_phase2,b_phase2,r,initial_cost,bv);
    final_cost = C'*X_opt;
    
    
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

    diff1 = f(:) - matrix_X_opt*onesMatrixf;
    diff2 = g(:) - (onesMatrixg'*matrix_X_opt)';
    maxDiff1 = max(diff1);
    maxDiff2 = max(diff2);
    
end

