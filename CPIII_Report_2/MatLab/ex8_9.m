function ex8_9(n)
% ex8-9 : switch between the differents exercice of the report 2.
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

% Session 8-9
%% Ex 1
function Exercice_1()

fprintf('Exercice 1 :\n');

% power
fprintf('\npower\n');
test_power

% ipower
fprintf('\nipower\n');
test_ipower

% rq
fprintf('\nrq\n');
test_rq

end

%% Ex 2
function Exercice_2()
fprintf('Exercice 2 :\n');

% 2.1
% fprintf('2.1\n');

% Define a new function for an ideal Hamiltonian
h_pure = @(N) diag(ones(1,N-1),1) + diag(ones(1,N-1),-1);
% Plot the Hamiltonian
imagesc(h_pure(10));

% 2.2
fprintf('2.2\n');

N = 10:5:20;

for i = 1:length(N)
    H = h_pure(N(i));
    figure
    imagesc(H);
    [vec, ~] = eig(H);
    figure
    imagesc(abs(vec.^2));
    
    % 2.3
    fprintf('2.3\n');
    [vec_min, ~] = eigs(H, 1, 'smallestreal');
    vec_min = conj(vec_min).*vec_min;
    
    H(1,N(i)) = 1;
    H(N(i),1) = 1;
    figure
    imagesc(H);
    [vec_min_2, ~] = eigs(H, 1, 'smallestreal');
    vec_min_2 = conj(vec_min_2).*vec_min_2;
    figure
    plot(1:N(i), vec_min, 1:N(i), vec_min_2)
    legend('Without hopping term', 'With hopping term')
    xlabel('x')
    ylabel('Probability')
    
    
    % 2.4
    fprintf('2.4\n');
    
    % Define a disorder to the Hamiltonian
    h_disorder = @(N,amplitude) diag(amplitude*round(rand(N,1))-amplitude/2);
    H = H + h_disorder(N(i), 1);
    figure
    imagesc(H);
    
    
    % 2.5
    fprintf('2.5\n');
    for j = 0:0.5:2
        H = h_pure(N(i)) + h_disorder(N(i),j);
        H(1,N(i)) = 1;
        H(N(i),1) = 1;
        [vec, ~] = eig(H);
        figure
        imagesc(abs(vec.^2));
    end
end

% 2.6
fprintf('2.6\n');

N = 100;
H = h_pure(N);
H(1,N) = 1;
H(N,1) = 1;
H_x = H + h_disorder(N, 1);
[vec_min_2, ~] = eigs(H, 1, 'smallestreal');
[vec_min_x, ~] = eigs(H_x, 1, 'smallestreal');
disp(sigma(vec_min_2));
disp(sigma(vec_min_x));
vec_min_2 = conj(vec_min_2).*vec_min_2;
vec_min_x = conj(vec_min_x).*vec_min_x;
figure
plot(1:N, vec_min_2, 1:N, vec_min_x)
legend('Without disorder', 'With disorder')
xlabel('x')
ylabel('Probability')

% 2.7
fprintf('2.7\n');

a = 0:10;
incert = zeros(length(a),1);
for i = 1:length(a)
    H = h_pure(N) + h_disorder(N,a(i));
    H(1,N) = 1;
    H(N,1) = 1;
    [vec, ~] = eig(H);
    for j = 1:length(vec)
        incert(i) = incert(i) + sigma(vec(:,j))/length(vec);
    end
end
figure
plot(a, incert)
xlabel('disorder')
ylabel('uncertainty')

% 2.8
fprintf('2.8\n');

a = 0.6;
N = 5:15:800;
incert = zeros(length(N),1);
for i = 1:length(N)
    H = h_pure(N(i)) + h_disorder(N(i),a);
    H(1,N(i)) = 1;
    H(N(i),1) = 1;
    [vec, ~] = eig(H);
    for j = 1:length(vec)
        incert(i) = incert(i) + sigma(vec(:,j))/length(vec);
    end
end
figure
plot(N, incert)
xlabel('size of matrix')
ylabel('uncertainty')
    
end

function s = sigma(psi)
N = length(psi);
% Define a coordinate operator
% N ? size of the Hamiltoian = number of atoms
x = diag(1:N);

s = sqrt(conj(psi)' * x.^2 * psi - (conj(psi)' * x * psi)^2);
end

%% Ex 3
function Exercice_3()
fprintf('Exercice 3 :\n');

% 3.1
fprintf('3.1\n');
test_j

% 3.2
fprintf('\n3.2\n');

s = 25;
precision = 5;
time_j  = zeros(s, 1);
time_cj = zeros(s, 1);
rotations_j  = zeros(s, 1);
rotations_cj = zeros(s, 1);
for i = 1:s
    for j = 1:precision
        clearvars m r t;
        m = rmg(1:i);
    
        tic;
        [~, r] = eig_j(m);
        t = toc;
        time_j(i) = time_j(i) + t;
        rotations_j(i) = rotations_j(i) + r;
    
        tic;
        [~, r] = eig_cj(m);
        t = toc;
        time_cj(i) = time_cj(i) + t;
        rotations_cj(i) = rotations_cj(i) + r;
    end
end

figure
plot(1:s, time_j./precision, 1:s, time_cj./precision)
legend('Time for j', 'Time for cj')
xlabel('size of matrix')
ylabel('Time for Jacobi calculation [s]')

figure
plot(1:s, rotations_j, 1:s, rotations_cj)
legend('Nb rotations j', 'Nb rotations cj')
xlabel('size of matrix')
ylabel('Number of rotations for Jacobi calculation')

end

%% Ex 4
function Exercice_4()
fprintf('Exercice 4 :\n');
% 4.1
fprintf('4.1\n');

V0 = -0.5;
r  = 5e-9;
eV = 1.60217656535e-19;
a  = 5*r;
a_tab = (1:0.01:9) .* r;
% N  = 30;
N  = 100;
% N  = 1000;
N_tab = [30 100 1000];
% N_tab = 30:10:160;
t = zeros(2, length(N_tab));
% b = zeros(2, length(N_tab), 'uint32');

for i = 1:length(N_tab)
    try
        [H_s, t_s] = H_sparse(N_tab(i), a, r);
%         disp(t_s);
        t(1, i) = t_s;
%         info = whos(H_s);
%         whos H_s;
%         b(1, i) = uint32(info.bytes);
        
        [H_f, t_f] = H_full(N_tab(i), a, r);
%         disp(t_f);
        t(2, i) = t_f;
%         info = whos(H_f);
%         whos H_f;
%         b(2, i) = uint32(info.bytes);
        
        disp(sparse(H_f-H_s));
    catch ME
        t(2, i) = 0;
%         b(2, i) = 0;
        disp(ME);
    end
end

figure
plot(N_tab, t(1,:), N_tab, t(2,:))
legend('Time for sparse', 'Time for full')
xlabel('size of matrix')
ylabel('Time for Hamiltonian creation [s]')

% figure
% plot(N_tab, b(1,:), N_tab, b(2,:))
% legend('Bytes for sparse', 'Bytes for full')
% xlabel('size of matrix')
% ylabel('Bytes for Hamiltonian creation')


% 4.2
fprintf('4.2\n');

[vec, val] = eigs(H_sparse(N, a, r), N, 'smallestreal');
disp(val(1:3,1:3));

counter = 0;
for i = 1:N
    if val(i,i) < 0
        counter = counter + 1;
    end
end
disp(counter);

phi_0 = zeros(N);
phi_1 = zeros(N);
phi_2 = zeros(N);

for i = 1:N^2
    phi_0(i) = vec(i, 1);
    phi_1(i) = vec(i, 2);
    phi_2(i) = vec(i, 3);
end

figure
pcolor(phi_0);
figure
pcolor(phi_1);
figure
pcolor(phi_2);

proba_0 = conj(phi_0).*phi_0;
proba_1 = conj(phi_1).*phi_1;
proba_2 = conj(phi_2).*phi_2;

figure
pcolor(proba_0);
figure
pcolor(proba_1);
figure
pcolor(proba_2);

V = reshape(V_pot(N, a, r)./V0/eV, [N,N]);

proba_0 = sum(sum(proba_0 .* V));
proba_1 = sum(sum(proba_1 .* V));
proba_2 = sum(sum(proba_2 .* V));

disp(proba_0);
disp(proba_1);
disp(proba_2);


% 4.3
fprintf('4.3\n');

[vec, ~] = eigs(H_sparse(N, a, r), 1, 1);

vec = reshape(vec, [N, N]);

figure
pcolor(vec);

proba_0 = conj(vec).*vec;

figure
pcolor(proba_0);

proba_0 = sum(sum(proba_0 .* V));

disp(proba_0);

prob = zeros(1,length(a_tab));

for i = 1:length(a_tab)
    [phi, val] = eigs(H_sparse(N, a_tab(i), r), 1, 1);
%     disp(val);
    
    phi = reshape(phi, [N,N]);
    
    proba = conj(phi).*phi;
    
    V = reshape(V_pot(N, a_tab(i), r)./V0/eV, [N,N]);
    probb = sum(sum(proba .* V));
    
    prob(i) = probb;
end

figure
plot(a_tab ./ r, prob)
xlabel('a/r')
ylabel('Probabitity of finding particle inside the quantum well')

end

function [H, t] = H_sparse(N, a, r)
% sparse matrix
eV = 1.60217656535e-19;
hb = 1.05457172647e-34;
me = 9.10938291e-31;
dr = a/(N-1);

tic;
matrix_sparse = sparse(1:N, 1:N, -4);
matrix_sparse = matrix_sparse + sparse(2:N, 1:N-1, 1, N, N) + sparse(1:N-1, 2:N, 1, N, N);

H = spalloc(N^2, N^2, 5*N^2);
for i = 1:N
    H((i-1)*N+1:i*N, (i-1)*N+1:i*N) = matrix_sparse;
end
H = H + sparse(N+1:N^2, 1:N^2-N, 1, N^2, N^2) + sparse(1:N^2-N, N+1:N^2, 1, N^2, N^2);
H = H ./ dr^2;

H = -hb^2 / (2 * me) * H + sparse(1:N^2, 1:N^2, V_pot(N, a, r));
H = H ./ eV;
t = toc;
end

function [H, t] = H_full(N, a, r)
% standard matrix
eV = 1.60217656535e-19;
hb = 1.05457172647e-34;
me = 9.10938291e-31;
dr = a/(N-1);

tic;
matrix = -4 * eye(N);
matrix = matrix + diag(ones(N-1,1),1) + diag(ones(N-1,1), -1);

H = zeros(N^2);
for i = 1:N
    H((i-1)*N+1:i*N, (i-1)*N+1:i*N) = matrix;
end
H = H + diag(ones(N^2-N,1), N) + diag(ones(N^2-N,1), -N);
H = H ./ dr^2;

H = -hb^2 / (2 * me) * H + diag(V_pot(N, a, r));
H = H ./ eV;
t = toc;
end

function V = V_pot(N, a, r)
% V for sparse & full
eV = 1.60217656535e-19;
V0 = -0.5;
dr = a/(N-1);
Nbr = r/dr;

V = zeros(N^2, 1);
for i = 1:N^2
    x = (N-1)/2 - mod(i-1, N);
    y = (N-1)/2 - floor((i-1)/N);
    if x^2 + y^2 < Nbr^2
        V(i) = V0 * eV;
    end
end
end
