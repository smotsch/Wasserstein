%Author: Sebastien Motsch (smotsch@asu.edu).

%This code takes takes a tableau and reduces it by Gaussian eleminitation
%until it arrives at an optimal solution.

function [A_up,b_up,X,bv]=MinimizeTableau(A,b,r,InitialCost,bvInit)
% Solve the minimization problem:
%     argmin C⋅X
%       with AX = b, X≥0
% The algorithm needs to start from a basic solution:
%       A = ['Id' | A_b^-1⋅A_n]
% Input:
%     . A : m×n matrix, m constraints, n free variables
%     . b : m vector, value of the constraint
%     . r : gradient of the cost in the non-basic variable
%             r = ∇_x_N (C⋅X) =c_n -  (A_b^-1 ⋅ A_n)^T c_b
% Output:
%     . Tableau: update tableau after row manipulation
%     . X  : optimal solution (argmin)
%     . bv : indices of the basic solutions
%     
    
% init
    [m n]= size(A);
    Tableau = [A b;
               r -InitialCost];
     
    if exist('bvInit', 'var')
        bv = bvInit;
    else
        bv = zeros(m,1);
    end
    %-------------------------%
    %   loop  minimization    %
    %-------------------------%
    tic
    for iloop = 1:n
        [valC j0] = min(Tableau(m+1,1:n)); % j0=ebv 'entering basic variable'
        if valC < -1d-10
            listRatio = Tableau(1:m,end)./Tableau(1:m,j0);
            listRatio = listRatio + (Tableau(1:m,j0)<=0)*1d10;
            [valRatio i0] = min(listRatio);                  % i0 lbv 'leaving basic variable'
            %% update basic variables
            bv(i0) = j0;
        else
            %fprintf('max!\n')
            break;
        end
        Tableau(i0,:) = Tableau(i0,:)/Tableau(i0,j0);     % Tableau(i0,j0) = 1
        for i = 1:(m+1)
            if (i ~= i0)
                Tableau(i,:) = Tableau(i,:) - Tableau(i,j0)*Tableau(i0,:);
            end
        end
        %Tableau(end,:);
    end
    timeAlgo=toc;
    %-------------------------%
    %-------------------------%
    % final

    jNotZero = (bv>0);
    X = zeros(n,1);
    X(bv(jNotZero)) = Tableau(jNotZero,end);
    % update A,b
    A_up = Tableau(1:m,1:n);
    b_up = Tableau(1:m,end);
    % print
    fprintf('Time algo: %f\n',timeAlgo)
    fprintf('Cost (before/after): %f / %f\n',InitialCost,-Tableau(end,end))
    