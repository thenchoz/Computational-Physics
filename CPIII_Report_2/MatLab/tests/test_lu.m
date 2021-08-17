% This script tests lu_decomposition.m

fprintf('Test 1: Identity matrix ...');
for n=1:1:10
    A=eye(n);
  [l u p]=lu_decomposition(A);
delta=p*A-l*u;
assert(trace(delta'*delta)<1e-10 && trace(triu(l,1)'*triu(l,1)) <1e-10 && trace(tril(u,-1)'*tril(u,-1)) <1e-10 );

end
fprintf('\tpassed\n');

fprintf('Test2 ...')
A=[1 2 3 4; 5 6 7 8; 6 8 10 9; 11 12 15 14];
[l u p]=lu_decomposition(A);
delta=p*A-l*u;
assert(trace(delta'*delta)<1e-10 && trace(triu(l,1)'*triu(l,1)) <1e-10 && trace(tril(u,-1)'*tril(u,-1)) <1e-10 );
fprintf('\tpassed\n');

fprintf('Test 3: Random real matrices ...')
for n=1:1:10
A=rand(n);
[l u p]=lu_decomposition(A);
delta=p*A-l*u;
assert(trace(delta'*delta)<1e-10 && trace(triu(l,1)'*triu(l,1)) <1e-10 && trace(tril(u,-1)'*tril(u,-1)) <1e-10 );
end
fprintf('\tpassed\n');

fprintf('Test 4: Random complex matrices ...')
for n=1:1:10
A=rand(n)+1i*rand(n);
[l u p]=lu_decomposition(A);
delta=p*A-l*u;
assert(trace(delta'*delta)<1e-10 && trace(triu(l,1)'*triu(l,1)) <1e-10 && trace(tril(u,-1)'*tril(u,-1)) <1e-10 );
end
fprintf('\tpassed\n')

fprintf('All tests passed!\n');