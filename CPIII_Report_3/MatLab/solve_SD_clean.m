function X = solve_SD(A, b)
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

X = rand(length(A), 1);
r = b - A*X;

epsilon = 1e-12;

while norm(r) > epsilon
    alpha = r'*r / (r' * A * r);
    X = X + alpha * r;
    r = b - A*X;
end
end