%This script tests solve_CG.m
%uses rmg.m

fprintf('Test 1: Identity matrix ...');
for n=1:1:10
    A=eye(n);
    X=rand(n,1);
    B=A*X;

assert(all(abs(solve_CG(A,B) - X)<1e-8));
end
fprintf('\tpassed\n');

fprintf('Test 2: Random matrices')
for n=1:1:5
   A=rand(n);
   A=A'*A;
   X=rand(n,1);
   B=A*X;
   
   assert(norm(solve_CG(A,B)-X) <1e-8);
      
end

fprintf('\tpassed\n');

fprintf('Test 3: Random matrices')
for n=1:1:5
    A=rmg(rand(n));
    X=n*rand(n,1);
    B=A*X;
    
    
assert(norm(solve_CG(A,B)-X)<1e-8);
end

fprintf('\tpassed\n');


fprintf('All tests passed!\n');