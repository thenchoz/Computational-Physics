function [ L, U, P ] = lu_decomposition( A )
% lu_decomposition computes the LU decomposition with pivoting 
% for a square matrix A
% Return lower L, upper U and permutation matrix P such that P*A=L*U  

% Initializing variables
N = length(A);
L = eye(N,N);
U = A;
P = eye(N,N);

for k = 1:N-1
    % Checking for the max value on the column
    [~, r] = max(abs(U(k:N,k)));
    r = r + k - 1;
    
    if r ~= k
        % if not the current row, permutation
        temp = U(k,:);
        U(k,:) = U(r,:);
        U(r,:) = temp;
        
        temp = P(k,:);
        P(k,:) = P(r,:);
        P(r,:) = temp;
        
        temp = L(k,1:k-1);
        L(k,1:k-1) = L(r,1:k-1);
        L(r,1:k-1) = temp;
    end
    
    for i = k+1:N
        L(i,k) = U(i,k)/U(k,k);
        for j = k:N
            U(i,j) = U(i,j) - L(i,k) * U(k,j);
        end
    end
end
end

