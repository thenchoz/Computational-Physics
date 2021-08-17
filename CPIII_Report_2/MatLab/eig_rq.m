function [ vec , val ] = eig_rq ( input_matrix , target )
% eig_rq:
%     computes the closest eigenvector and eigenvalue of a given matrix .
%
% Arguments:
% 
%       input_matrix (2D complex Hermitian matrix):
%                                    matrix for the eigenvalue problem ;
%       target (real scalar): an estimation to the eigenvalue;
%
% Returns:
%       a right eigenvector and the corresponding eigenvalue of a matrix.

if input_matrix == diag(diag(input_matrix))
    % case input_matrix is diagonal, avoiding singular matrix
    [~, i] = min(abs(diag(input_matrix) - target));
    val = input_matrix(i,i);
    vec = zeros(length(input_matrix), 1);
    vec(i) = 1;
else
    % initialisation
    id = eye(size(input_matrix));
    
    vec = rand(length(input_matrix),1);
    vec = vec/norm(vec);
    val = target;
    
    % to enter the loop
    temp = target + 1;
    
    epsilon = 1e-12;
    
    while abs(val - temp) > epsilon
        temp = val;
        vec = (input_matrix - temp * id) \ vec;
        vec = vec/norm(vec);
        val = vec'*input_matrix*vec;
    end
end
end