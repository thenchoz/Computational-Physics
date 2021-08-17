function [ x ] = solve_ls( A, b, lu_choice)
% solve_ls solves the system of linear equations A*x=b, for a square matrix
% A and a column vector b, using the LU decomposition and backward 
% substitution


if(nargin == 2)
    lu_choice = 0;
end

switch lu_choice
    case 1
        [L U] = lu_decomposition_withoutP(A);
    case 2
        [L U P] = lu(A);
        % permuting b because LU = PA
        b = P*b;
    otherwise
        [L U P] = lu_decomposition(A);
        % permuting b because LU = PA
        b = P*b;
end

N = length(A);

% resolving Ly = b
% usualy, L(i,i) = 1 but...
y = b;
y(1) = 1 / L(1,1) * b(1);
for i = 2:N
    y(i) = 1 / L(i,i) * (b(i) - sum(L(i,1:i-1) * y(1:i-1)));
end

% resolving Ux = y
x = y;
x(N) = 1 / U(N,N) * y(N);
for j = 1:N-1
    i = N-j;
    x(i) = 1 / U(i,i) * (y(i) - sum(U(i,i+1:N) * x(i+1:N)));
end
end

