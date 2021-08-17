function [val, rotations] = eig_cj(input_matrix)
% eig_cj computes the eigenvalues of a matrix.
% 
% Arguments:
%     input_matrix (2D real symmetric matrix):
%                   matrix for the eigenvalue problem;
%     
% Returns : an array with eigenvalues.

A = input_matrix;
N = length(A);
% *N^2 to have all value close to epsilon for big matrix
epsilon = 1e-12 * N^2;
rotations = 0;

while off(A) > epsilon
    for p = 1:N-1
        for q = p+1:N
            rotations = rotations + 1;
            [c, s] = cs_find(A, p, q);
            
            % J'AJ without matrix multiplication
            temp = A(:, p);
            A(:, p) = c * temp - s * A(:,q);
            A(:, q) = s * temp + c * A(:,q);
            
            temp = A(p, :)';
            A(p, :) = c * temp' - s * A(q,:);
            A(q, :) = s * temp' + c * A(q,:);
        end
    end
end

val = diag(A);

end

function o = off(A)
A = A.^2;
o = sqrt(sum(sum(A)) - trace(A));
end

function [c, s] = cs_find(A, p, q)
if A(p, q) ~= 0
    tau = (A(q, q) - A(p, p)) / (2 * A(p, q));
    if tau >= 0
        t = -tau + sqrt(1 + tau^2);
    else
        t = -tau - sqrt(1 + tau^2);
    end
    c = 1 / sqrt(1 + t^2);
    s = t * c;
else
    c = 1;
    s = 0;
end
end