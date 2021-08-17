function X = solve_CG(A, b)
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

X = rand(length(A), 1);
r = b - A*X;
d = r;

epsilon = 1e-12;

while norm(r) > epsilon
    alpha = r'*r / (d' * A * d);
    X = X + alpha * d;
    r_new = r - alpha * A * d;
    alpha = r_new'*r_new / (r'*r);
    r = r_new;
    d = r + alpha * d;
end
end