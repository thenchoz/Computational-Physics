function [L, U] = lu_decomposition_withoutP( A )

% lu_decomposition computes the LU decomposition without pivoting 
% for a square matrix A
% Return lower L and upper U such that A=L*U  

N = length(A);
L = eye(N,N);
U = A;

for k = 1:N-1
    for i = k+1:N
        L(i,k) = U(i,k)/U(k,k);
        for j = k:N
            U(i,j) = U(i,j) - L(i,k) * U(k,j);
        end
    end
end

end

