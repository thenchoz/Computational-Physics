%This script tests solve_ls.m

fprintf('Test 1: Identity matrix ...');
for n=1:1:10
    A=eye(n);
    X=rand(n,1);
    B=A*X;

assert(all(abs(solve_ls(A,B) - X)<1e-12));
end
fprintf('\tpassed\n');

fprintf('Test 2: Random real matrices')
for n=1:1:10
   A=rand(n);
   X=rand(n,1);
   B=A*X;
   
   assert(all(abs(solve_ls(A,B)-X)<1e-12));
      
end
fprintf('\tpassed\n');
   
fprintf('Test 3: Random imaginary matrices')
for n=1:1:10
    
    A=rand(n)+1i*rand(n);
    X=rand(n,1)+1i*rand(n,1);
    B=A*X;
    
    assert(all(abs(solve_ls(A,B)-X)<1e-12));
end
fprintf('\tpassed\n');

fprintf('All tests passed!\n');