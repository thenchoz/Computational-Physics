function ex12_13(n)
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
    case 3
        Exercice_3()
end
end

% Session 12-13
function Exercice_1()
fprintf('Exercice 1 :\n');
% 1.1
fprintf('1.1\n');

A = [
    1 2;
    3 4;
    5 6];

b = [1;1;2];

x = mypinv(A,b);
disp(x);

% 1.2
fprintf('1.2\n');

N = 100;
vicinity(1,:) = linspace(x(1) - 1, x(1) + 1, N);
vicinity(2,:) = linspace(x(2) - 1, x(2) + 1, N);

test = zeros(N);
for i = 1:N
    for j = 1:N
        p(1) = vicinity(1,i);
        p(2) = vicinity(2,j);
        test(i,j) = norm(A*p' - b);
    end
end

% figure
% surf(vicinity(1,:), vicinity(2,:), test);

figure
hold on
contourf(vicinity(1,:), vicinity(2,:), test, 'LevelList', -20:0.1:20);
plot(x(1), x(2), 'r*');
hold off

end


function Exercice_2()
fprintf('Exercice 2 :\n');
% 2.1
fprintf('2.1\n');

C1 = [
    0.3718 0.3464 0.4041;
    0.3646 0.0156 0.3409;
    0.1925 0.1359 0.1556;
    0.3216 0.0227 0.4663;
    0.0840 0.0477 0.0169];

C2 = [
    0.2276 0.4382 0.1554;
    0.0870 0.1676 0.0594;
    0.2921 0.5625 0.1995;
    0.0992 0.1910 0.0677;
    0.1967 0.3789 0.1344];

[U1, S1, V1] = svd(C1);
[U2, S2, V2] = svd(C2);

disp(S1);
disp(S2);

% 2.2
fprintf('2.2\n');

% Psi_A = U1 * S1;
% Psi_B = S1 * V1;
% 
% disp(Psi_A);
% disp(Psi_B);

Psi_A = U2(:,1);
Psi_B = V2(:,1);

disp(Psi_A);
disp(Psi_B);


end


function Exercice_3()
fprintf('Exercice 3 :\n');
% 3.1
fprintf('3.1\n');

data_Mpc   = [0.032 0.034 0.214 0.263 0.275 0.275 0.45 0.5 0.5 0.63 0.8 0.9 0.9 0.9 0.9 1.0 1.1 1.1 1.4 1.7 2.0 2.0 2.0 2.0];
data_speed = [170 290 -130 -70 -185 -220 200 290 270 200 300 -30 650 150 500 920 450 500 500 960 500 850 800 1090];

Matrix(:,1) = data_Mpc;
Matrix(:,2) = ones(size(data_Mpc));
ab = mypinv(Matrix, data_speed');
Matrix(:,2) = zeros(size(data_Mpc));
a_wb = mypinv(Matrix, data_speed');

fit = ab(1) * data_Mpc + ab(2);
fit_wb = a_wb(1) * data_Mpc;

disp(ab(1));
disp(a_wb(1));


% 3.2
fprintf('3.2\n');

figure
hold on
plot(data_Mpc, fit)
plot(data_Mpc, fit_wb)
plot(data_Mpc, data_speed, '.')
xlabel('Distance [Mpc]')
ylabel('Speed [km/h]')
hold off

end

function X = mypinv(A, b)

[U, S, V] = svd(A);
d = U'*b;
if S(2,2) == 0
    z = [1/S(1,1) * d(1); 0];
else
    z = [1/S(1,1) * d(1); 1/S(2,2) * d(2)];
end
X = V * z;

end