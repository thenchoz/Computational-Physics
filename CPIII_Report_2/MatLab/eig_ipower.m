function [ vec , val ] = eig_ipower( input_matrix , target )
% eig_ipower:
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

% initialisation
id = eye(size(input_matrix));

vec = rand(length(input_matrix),1);
vec = vec/norm(vec);
val = vec'*input_matrix*vec;

% to enter the while
temp = val + 1;

epsilon = 1e-15;

while abs(val - temp) > epsilon
    temp = val;
    vec = (input_matrix - target * id) \ vec;
    vec = vec/norm(vec);
    val = vec'*input_matrix*vec;
end
end