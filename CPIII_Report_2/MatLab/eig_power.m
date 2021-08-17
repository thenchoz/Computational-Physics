function [ vec , val ] = eig_power(input_matrix)
% eig_power: computes an eigenvector and eigenvalue of a given matrix .
%
% Arguments:
%       input_matrix (2D complex Hermitian matrix):
%                                    matrix for the eigenvalue problem ;
%
% Returns:
%       a right eigenvector and the corresponding eigenvalue of a matrix.

% initialisation
vec = rand(length(input_matrix),1);
vec = vec/norm(vec);
val = vec'*input_matrix*vec;

% to enter the while
temp = val + 1;

epsilon = 1e-15;

while abs(val - temp) > epsilon
    temp = val;
    vec = input_matrix * vec;
    vec = vec/norm(vec);
    val = vec'*input_matrix*vec;
end
end