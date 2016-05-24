% This code computes the distance between sets of points in R^2, R^3 by calling the
% Tableau code.
% Authors: Shane Lubold (shane.Lubold@asu.edu) & Sebastien Motsch (smotsch@asu.edu)

%Solve the optimal transport problem between two sets of points in R^2:
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
clear; close all; clc;

response = input('Do you want to optimize in 2D or 3D? ','s');
switch response
    case '2D'
        
    % A) initial configuration
    %--------------------------------

        Xpos = randn(4,2); % Starting points
        Ypos = randn(4,2) + ones(4,1)*[3 0]; % Ending points
         Xpos = [0 2 3; 1 2 5]';
         Ypos = [0 1 5; 2 3 3]';
%         
        nX = length(Xpos(:,1));
        nY = length(Ypos(:,1));

    case '3D'

        Xpos = randn(5,3); % Starting points
        Ypos = randn(5,3) + ones(5,1)*[5 0 1]; % Ending points
        nX = length(Xpos(:,1));
        nY = length(Ypos(:,1));
             
end
        

% B) constraint:
%---------------
% We are lookin for a transition matrix  '(x_ij)' with
%   x_ij piece of i going to j
% '(x_ij)' is a nX×nY matrix.  
% Thus,
% ∑_i x_ij = 1/nX   -> each Xpos particle 'weight' 1/nX
% ∑_j x_ij = 1/nY   -> each Ypos particle 'weight' 1/nY
% We transform the matrix '(a_ij)' into a vector 'X'
b = [1/nX*ones(nX,1);
     1/nY*ones(nY,1)];
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

matrix_C = zeros(nX,nY);

switch response
        case '2D'
for i=1:nX
    for j=1:nY
        %matrix_C(i,j) = (Ypos(j,1)-Xpos(i,1))^2 + (Ypos(j,2)-Xpos(i,2))^2;
        matrix_C(i,j) = sqrt( (Ypos(j,1)-Xpos(i,1))^2 + (Ypos(j,2)-Xpos(i,2))^2 );
    end
end
    C = matrix_C(:);

        case '3D'
for i=1:nX
    for j=1:nY
        %matrix_C(i,j) = (Ypos(j,1)-Xpos(i,1))^2 + (Ypos(j,2)-Xpos(i,2))^2;
        matrix_C(i,j) = sqrt( (Ypos(j,1)-Xpos(i,1))^2 + (Ypos(j,2)-Xpos(i,2))^2 + (Ypos(j,3)-Xpos(i,3))^2 );
    end
end
    C = matrix_C(:);        
end
           
%---------------------------------------------------------------%
%-----------------         Resolution            ---------------%
%---------------------------------------------------------------%

%----- Phase 1 -----%
%-------------------%
fprintf('---------------------------------------\n')
fprintf('Phase 1: find a basic feasible solution\n')
fprintf('---------------------------------------\n')
if (nX==nY)
    % simple basic feasible solution: 
    %    Xi is sent to Yi
    %   and Xn sent to no-one else
    %  -> better to use the Hungarian algorithm in this case
    %     http://cstheory.stackexchange.com/questions/2707/justification-for-the-hungarian-method-kuhn-munkres
    tic
    matrix_X = diag(ones(nX,1));
    X        = 1/nX*matrix_X(:);
    A_phase2 = A;
    b_phase2 = b;
    for i=1:(nX-1)
        A_phase2(i+nX,:) = A_phase2(i+nX,:) - A_phase2(i,:);
        b_phase2(i+nX,:) = b_phase2(i+nX,:) - b_phase2(i,:);
    end
    for i=1:(nX-1)
        A_phase2(nX,:) = A_phase2(nX,:) - A_phase2(nX+i,:);
        b_phase2(nX,:) = b_phase2(nX,:) - b_phase2(nX+i,:);
    end    
    bv         = zeros(2*nX-1,1);
    bv(1:nX)   = [1 1+(nX+1)*(1:(nX-1))];
    bv((nX+1):end) = nX^2-[(nX-1):(-1):1];
    timePhase1=toc
else
    [A_phase2,b_phase2,X,bv] = MinimizeTableau(A,b,-sum(A),-sum(b));
end
% reduce gradient 'r'
r = C';
for i=1:length(bv)
    % find associated basic variable
    j = bv(i);
    if (j>0)
        % update cost
        %  ∇_x_N r = c_N' - c_b'⋅(A_b^-1 A_n)
        %          = c'   - c_b'⋅(A_b^-1 A)
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
switch response
    case '2D'
% plot the result
matrix_X_opt = reshape(X_opt,nY,nX)';
clf
hold on

plot(Xpos(:,1),Xpos(:,2),'kx',Ypos(:,1),Ypos(:,2),'ko')
grid on
legend('Starting Points','Ending Points');

for i=1:nX
    for j=1:nY
        if matrix_X_opt(i,j)>1d-10
            quiver([Xpos(i,1),0],[Xpos(i,2),0], ...
                   [Ypos(j,1)-Xpos(i,1),0],[Ypos(j,2)-Xpos(i,2),0],'k--','autoScale','off')
        end
    end
end

title('Optimal Solution')
hold off

xlabel('x');
ylabel('y');

    case '3D'
        matrix_X_opt = reshape(X_opt,nY,nX)';
clf
hold on
plot3(Xpos(:,1),Xpos(:,2),Xpos(:,3),'bo',Ypos(:,1),Ypos(:,2),Ypos(:,3),'rx')
grid on
legend('Starting Points','Ending Points');

for i=1:nX
    for j=1:nY
        if matrix_X_opt(i,j)>1d-10
            quiver3([Xpos(i,1),0],[Xpos(i,2),0],[Xpos(i,3),0], ...
                   [Ypos(j,1)-Xpos(i,1),0],[Ypos(j,2)-Xpos(i,2),0],[Ypos(j,3)-Xpos(i,3),0],'k--','autoScale','off')
        end
    end
end

title('Optimal Solution')
hold off

xlabel('x');
ylabel('y');
zlabel('z');
view([50,-50,50]);
end
