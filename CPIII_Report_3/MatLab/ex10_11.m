function ex10_11(n)
% ex10-11 : switch between the differents exercice of the report 3.
%
% Arguments:
%       input interger: number of exercice;
%
% Returns : nothing.

switch n
    case 1
        Exercice_1()
    case 2
        Exercice_2()
end
end

% Session 10-11
function Exercice_1()

fprintf('Exercice 1 :\n');
% 1.1
fprintf('1.1\n');

test_solve_SD;
test_solve_CG;

% 1.2
fprintf('1.2\n');

matrix = matfile('Matrices.mat');
[~, i_SD, e_SD] = solve_SD(matrix.A1, matrix.b1);
[~, i_CG, e_CG] = solve_CG(matrix.A1, matrix.b1);

figure
plot(1:i_SD, e_SD, '.', 1:i_CG, e_CG, '.')
legend('SD method', 'CG method')
xlabel('iteration')
ylabel('error')

disp(cond(matrix.A1));

% 1.3
fprintf('1.3\n');

% [~, i_SD, e_SD] = solve_SD(matrix.A2, matrix.b2);
% [~, i_CG, e_CG] = solve_CG(matrix.A2, matrix.b2);
% 
% figure
% plot(1:i_SD, e_SD, '.', 1:i_CG, e_CG, '.')
% legend('SD method', 'CG method')
% xlabel('iteration')
% ylabel('error')

disp(cond(matrix.A2));

end


function Exercice_2()
fprintf('Exercice 2 :\n');
% 2.1
fprintf('2.1\n');

    function e = energy(X)
        N = length(X) / 3;
        X = reshape(X, 3, N);
        X = X';
        e = 0;
        
        for i = 1:N-1
            for j = i+1:N
                e = e + ULJ(norm(X(i,:) - X(j,:)));
            end
        end
    end

    function ulj = ULJ(x)
        sigma = 1;
        epsilon = 1;
        
        ulj = 4 * epsilon * ((sigma/x)^12 - (sigma/x)^6);
    end

step = 1e-8;

% 2.2
fprintf('2.2\n');

min_energy = 1;
thr = 1.13;

n = 4;
X = zeros(3 * n, 1);
% m2_x
X(4) = min_energy;
% m3_xy
X(7) = min_energy / 2;
X(8) = sqrt(3)/2 * min_energy;
% m4_xyz
X(10) = min_energy / 2;
X(11) = 1/sqrt(12) * min_energy;
X(12) = sqrt(5/12) * min_energy;

e_tetra = energy(X);
x_test = reshape(X, 3, n)';
visualize_molecule(x_test, thr, 1, 'tetrahedral init');
fprintf('tetrahedral init\n');
disp(e_tetra);

X_solve = solve_nLCG(@energy, step, X);
e_tetra = energy(X_solve);
X_solve = reshape(X_solve, 3, n)';
visualize_molecule(X_solve, thr, 2, 'tetrahedral final');
fprintf('tetrahedral final\n');
disp(e_tetra);



X = zeros(3 * n, 1);
% m2_x
X(4) = min_energy;
% m3_y
X(8) = min_energy;
% m4_xy
X(10) = min_energy;
X(11) = min_energy;

X_solve = solve_nLCG(@energy, step, X);
e_square = energy(X_solve);
X_solve = reshape(X_solve, 3, n)';
visualize_molecule(X_solve, thr, 3, 'square planar');
fprintf('square planar\n');
disp(e_square);




% 2.3
fprintf('2.3\n');

n = 5;
X = zeros(3 * n, 1);
% m2_x
X(4) = min_energy;
% m3_xy
X(7) = min_energy / 2;
X(8) = sqrt(3)/2 * min_energy;
% m4_xyz
X(10) = min_energy / 2;
X(11) = 1/sqrt(12) * min_energy;
X(12) = sqrt(5/12) * min_energy;
% m5_xyz
X(13) = min_energy / 2;
X(14) = 1/sqrt(12) * min_energy;
X(15) = -sqrt(5/12) * min_energy;

e_bipyramid = energy(X);
x_test = reshape(X, 3, n)';
visualize_molecule(x_test, thr, 4, 'trigonal bipyramid init');
fprintf('trigonal bipyramid init\n');
disp(e_bipyramid);

X_solve = solve_nLCG(@energy, step, X);
e_bipyramid = energy(X_solve);
X_solve = reshape(X_solve, 3, n)';
visualize_molecule(X_solve, thr, 5, 'trigonal bipyramid final');
fprintf('trigonal bipyramid final\n');
disp(e_bipyramid);



X = zeros(3 * n, 1);
% m2_x
X(4) = min_energy;
% m3_y
X(8) = min_energy;
% m4_xy
X(10) = min_energy;
X(11) = min_energy;
% m5_xyz
X(13) = min_energy / 2;
X(14) = min_energy / 2;
X(15) = 1/sqrt(2) * min_energy;

e_spyramid = energy(X);
x_test = reshape(X, 3, n)';
visualize_molecule(x_test, thr, 6, 'square pyramid init');
fprintf('square pyramid init\n');
disp(e_spyramid);

X_solve = solve_nLCG(@energy, step, X);
e_spyramid = energy(X_solve);
X_solve = reshape(X_solve, 3, n)';
visualize_molecule(X_solve, thr, 7, 'square pyramid final');
fprintf('square pyramid final\n');
disp(e_spyramid);




% 2.4
fprintf('2.4\n');

n = 6;
X = zeros(3 * n, 1);
% m1_x
X(1) = sqrt(2) * min_energy;
% m2_y
X(5) = sqrt(2) * min_energy;
% m3_x
X(7) = -sqrt(2) * min_energy;
% m4_y
X(11) = -sqrt(2) * min_energy;
% m5_z
X(15) = sqrt(2) * min_energy;
% m6_z
X(18) = -sqrt(2) * min_energy;

e_octahedral = energy(X);
x_test = reshape(X, 3, n)';
visualize_molecule(x_test, 3*thr, 8, 'octahedral init');
fprintf('octahedral init\n');
disp(e_octahedral);

X_solve = solve_nLCG(@energy, step, X);
e_octahedral = energy(X_solve);
X_solve = reshape(X_solve, 3, n)';
visualize_molecule(X_solve, 7*thr, 9, 'octahedral final');
fprintf('octahedral final\n');
disp(e_octahedral);



X = zeros(3 * n, 1);
% m1_x
X(1) = sin(pi/5) * min_energy;
% m2_xy
X(4) = cos(2/5*pi) * sin(pi/5) * min_energy;
X(5) = sin(2/5*pi) * sin(pi/5) * min_energy;
% m3_xy
X(7) = cos(2/5*pi) * sin(pi/5) * min_energy;
X(8) = -sin(2/5*pi) * sin(pi/5) * min_energy;
% m4_xy
X(10) = -cos(pi/5) * sin(pi/5) * min_energy;
X(11) = sin(pi/5) * sin(pi/5) * min_energy;
% m5_xy
X(13) = -cos(pi/5) * sin(pi/5) * min_energy;
X(14) = -sin(pi/5) * sin(pi/5) * min_energy;
% m6_z
X(18) = sqrt(1 - sin(pi/5)^2) * min_energy;

e_pentagonal = energy(X);
x_test = reshape(X, 3, n)';
visualize_molecule(x_test, 3*thr, 10, 'pentagonal pyramid init');
fprintf('pentagonal pyramid init\n');
disp(e_pentagonal);

X_solve = solve_nLCG(@energy, step, X);
e_pentagonal = energy(X_solve);
X_solve = reshape(X_solve, 3, n)';
visualize_molecule(X_solve, 3*thr, 11, 'pentagonal pyramid final');
fprintf('pentagonal pyramid final\n');
disp(e_pentagonal);


end