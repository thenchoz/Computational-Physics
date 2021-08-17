function ex6_7(n)
% ex6-7 : switch between the differents exercice of the report 2.
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
    case 3
        Exercice_3()
    case 4
        Exercice_4()
end
end

% Session 6-7
function Exercice_1()

fprintf('Exercice 1 :\n');
% 1.1
fprintf('1.1\n');

A = [1 2 3 4;-1 2 -3 4;0 1 -1 1;1 1 1 1];
B = [30 10 3 10];
X = A\B'

X = solve_ls(A, B', 0)

% 1.2

s = 50;
pres = 5;
diff = zeros(7, s);
for i = 1:s
    for j = 1:pres
        A = rand(i);
        tic;
        [l u p] = lu(A);
        t1 = toc;
        tic;
        [L U P] = lu_decomposition(A);
        t2 = toc;
        [K V] = lu_decomposition_withoutP(p * A);
        diff(1, i) = diff(1, i) + ((sum(sum(abs(l-L))) + sum(sum(abs(u-U)))) / i^2);
        diff(2, i) = diff(2, i) + ((sum(sum(abs(l-K))) + sum(sum(abs(u-V)))) / i^2);
        diff(3, i) = diff(3, i) + (sum(sum(abs(l*u - p*A))) / i^2);
        diff(4, i) = diff(4, i) + (sum(sum(abs(L*U - P*A))) / i^2);
        diff(5, i) = diff(5, i) + (sum(sum(abs(K*V - p*A))) / i^2);
        diff(6, i) = diff(6, i) + t1;
        diff(7, i) = diff(7, i) + t2;
    end
end

figure
plot(1:s, diff(1,:)./pres, 1:s, diff(3,:)./pres, 1:s, diff(5,:)./pres)
legend('Dif. between our L \& U and Matlab L \& U', 'Dif. between Matlab LU and PA', 'Dif. between our LU and PA')
xlabel('size of matrix')
ylabel('average error for each point')

figure
plot(1:s, diff(6,:)./pres, 1:s, diff(7,:)./pres)
legend('Matlab function', 'Our function')
xlabel('size of matrix')
ylabel('average time used to calculate [s]')

% 1.3
fprintf('1.3\n');

A = [1 2 0;2 4 8;3 -1 2];
[l u] = lu_decomposition_withoutP(A)
% division by 0 ?

% 1.4
fprintf('1.4\n');

[L U P] = lu_decomposition(A)
test_lu
fprintf('\n');

test_solve
fprintf('\n');

end


function Exercice_2()
fprintf('Exercice 2 :\n');
% 2.1
% fprintf('2.1\n');

syms I1 I2 I3 I4 I R1 R2 R3 R4 R V; % U1 U2 U3 U4 U;
% Ohm : Ui = Ri*Ii; U = R*I
% Noeuds : I1 = I3 + I; I2 + I = I4
% Mailles : U1 + U = U2; U3 = U + U4; V = U1 + U3 = U2 + U4

Matrix_Wheatstone = [
    1  0 -1   0  -1;
    0 -1  0   1  -1;
   -R1 R2 0   0  -R;
    0  0  R3 -R4 -R;
    R1 0  R3  0   0
    ];
Intensity = [I1; I2; I3; I4; I];
Result = [0; 0; 0; 0; V];
% Matrix_Wheastone * Intensity = Result

% 2.2
fprintf('2.2\n');

Sol = solve(Matrix_Wheatstone * Intensity == Result, I1, I2, I3, I4, I);
disp(Sol.I);

% 2.3
fprintf('2.3\n');

syms_R = [R1; R2; R3; R4; R];
value_R = [20000; 40000; 15000; 30000; 10000];
syms_V = V;
value_V  = 10;
Matrix_Wheatstone = subs(Matrix_Wheatstone, syms_R, value_R);
Result = subs(Result, syms_V, value_V);

disp(Matrix_Wheatstone);

Intensity_solved(:, 1) = solve_ls(Matrix_Wheatstone, Result, 0);
Intensity_solved(:, 2) = solve_ls(Matrix_Wheatstone, Result, 1);
Intensity_solved(:, 3) = solve_ls(Matrix_Wheatstone, Result, 2);
Intensity_sol_dif(:, 1) = abs(Intensity_solved(:,1) - Intensity_solved(:,3));
Intensity_sol_dif(:, 2) = abs(Intensity_solved(:,2) - Intensity_solved(:,3));

disp(Intensity_solved);
disp(Intensity_sol_dif);

Matrix_Wheatstone_2 = [
    1  0 -1   0  -1;
    0 -1  0   1  -1;
   -20 40 0   0  -10;
    0  0  15 -30 -10;
    20 0  15  0   0
    ];
Result_2 = [0; 0; 0; 0; 10];

disp(Matrix_Wheatstone_2);

Intensity_solved_2(:, 1) = solve_ls(Matrix_Wheatstone_2, Result_2, 0);
Intensity_solved_2(:, 2) = solve_ls(Matrix_Wheatstone_2, Result_2, 1);
Intensity_solved_2(:, 3) = solve_ls(Matrix_Wheatstone_2, Result_2, 2);
Intensity_sol_dif_2(:, 1) = abs(Intensity_solved_2(:,1) - Intensity_solved_2(:,3));
Intensity_sol_dif_2(:, 2) = abs(Intensity_solved_2(:,2) - Intensity_solved_2(:,3));

disp(Intensity_solved_2);
disp(Intensity_sol_dif_2);

end


function Exercice_3()
fprintf('Exercice 3 :\n');

% n1[FeS] + n2[NaBiO3] + n3[H2SO4] -> n4[Bi2(SO4)3] + n5[Fe2(SO4)3] + n6[Na2SO4] + n7[H2O]
syms Fe S Na Bi O H n1 n2 n3 n4 n5 n6 n7;
Element = [Fe; S; Na; Bi; O; H];
Coef = [n1 n2 n3 n4 n5 n6 n7];
Matrix_coeff_left = [
    n1 n1 0  0  0    0;
    0  0  n2 n2 3*n2 0;
    0  n3 0  0  4*n3 2*n3
    ];
Matrix_coeff_right = [
    0    3*n4 0    2*n4 12*n4 0;
    2*n5 3*n5 0    0    12*n5 0;
    0    n6   2*n6 0     4*n6 0;
    0    0    0    0    n7    2*n7
    ];
equ  = sum(Matrix_coeff_left) - sum(Matrix_coeff_right);
disp(equ*Element);

Sol = solve(equ == 0, n1, n2, n3, n4, n5, n6);%, n7);
multiple(1) = Sol.n1;
multiple(2) = Sol.n2;
multiple(3) = Sol.n3;
multiple(4) = Sol.n4;
multiple(5) = Sol.n5;
multiple(6) = Sol.n6;
multiple(7) = n7;

disp(Coef == multiple);

% ugly but no other idea yet
for i = 1:ceil(prod(n7./multiple))
    temp = subs(multiple, n7, i);
    if all(mod(temp, 1) == zeros(size(temp)))
        break
    end
end

disp(Coef == temp);

end


function Exercice_4()
fprintf('Exercice 4 :\n');

syms x1 x2 x3;
Matrix = [
    1     0 1;
    1.001 1 0;
    0    -1 1
    ];
Result1 = [2; 1; 1];
Result2 = [2.001; 1; 1];

result_x_1 = solve_ls(Matrix, Result1);
result_x_2 = solve_ls(Matrix, Result2);

disp(cond(Matrix));
disp(result_x_1);
disp(result_x_2);
disp(Matrix * result_x_1);
disp(Matrix * result_x_2);

end
