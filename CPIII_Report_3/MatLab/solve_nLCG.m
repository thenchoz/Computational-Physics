function [X] = solve_nLCG(f, step, x)
% solve_nLSD:
%     find the minimum potential position from the start position x and
%     function f
%
% Arguments:
%       f : function to minimise;
%       step : step for the gradient/hessian
%       x : init position
%
% Returns:
%       the n minima position in matrix X.

X = x;
r = - grad(f, step, X);
d = r;
i = 0;

epsilon = 1e-5 * length(x);

while norm(r) > epsilon
    i = i + 1;
    if mod(i, 50000) == 0
        fprintf(num2str(i));
        fprintf('\t:\t');
        fprintf(num2str(f(X)));
        fprintf('\t:\t');
        fprintf(num2str(norm(r)));
        fprintf('\n');
    end
%     alpha = 1;
%     while norm(alpha * d) > epsilon
        h = findH(f, step, x, d);
        if h > 0
            alpha = - grad(f, step, X)'*d / h;
            X = X + alpha * d;
        else
%             if f(X - 1e-3 * d ./ norm(d)) < f(X + 1e-3 * d ./ norm(d))
                X = X - 1e-3 * d ./ norm(d);
%             else
%                 X = X + 1e-3 * d ./ norm(d);
%             end
        end
%     end
    r_new = - grad(f, step, X);
    alpha = r_new'*r_new / (r'*r);
    r = r_new;
    d = r + alpha * d;
end
end

function g = grad(f, step, X)
    N = length(X);
    g = zeros(N, 1);
    D = zeros(N, 1);
    for i = 1:N
        D(i) = step / 2;
        g(i) = (f(X+D) - f(X-D)) / step;
        D(i) = 0;
    end
end

function [dhd] = findH(f, step, X, d)
% d matrix direction

dhd = (f(X + step*d) + f(X - step*d) - 2*f(X)) / step^2;
end

