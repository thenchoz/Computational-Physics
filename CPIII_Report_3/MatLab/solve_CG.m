function [X, i, e] = solve_CG(A, b)
% solve_SD:
%     solves the system of linear equations A*x=b, using the conjugate
%     gradient method
%
% Arguments:
%       A : square matrix;
%       b : column vector;
%
% Returns:
%       a solution vector x;
%       the number of iteration i;
%       all approximation X.

X = rand(length(A), 1);
r = b - A*X;
d = r;
i = 0;

% !! calcul of the real solution, only for returning the error
x_real = A\b;

epsilon = 1e-12;

while norm(r) > epsilon
    alpha = r'*r / (d' * A * d);
    X = X + alpha * d;
    r_new = r - alpha * A * d;
    alpha = r_new'*r_new / (r'*r);
    r = r_new;
    d = r + alpha * d;
    i = i + 1;
    e(i) = norm(X - x_real)/norm(x_real);
end
end