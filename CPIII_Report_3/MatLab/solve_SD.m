function [X, i, e] = solve_SD(A, b)
% solve_SD:
%     solves the system of linear equations A*x=b, using the steepest
%     descent method
%
% Arguments:
%       A : square matrix;
%       b : column vector;
%
% Returns:
%       a solution vector x;
%       the number of iteration i;
%       error e for each x.

X = rand(length(A), 1);
r = b - A*X;
i = 0;

% !! calcul of the real solution, only for returning the error
x_real = A\b;

epsilon = 1e-12;

if length(A) > 15
    % for big matrix, awoiding too long loop
    epsilon = epsilon * length(A);
end

while norm(r) > epsilon
    alpha = r'*r / (r' * A * r);
    X = X + alpha * r;
    r = b - A*X;
    i = i + 1;
    e(i) = norm(X - x_real)/norm(x_real);
end
end