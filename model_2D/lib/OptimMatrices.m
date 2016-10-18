function [A, b, C ] = OptimMatrices(f,g,x,y)

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


end

