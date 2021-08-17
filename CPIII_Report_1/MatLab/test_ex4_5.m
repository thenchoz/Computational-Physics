% This script tests myfft.m for correctness

fprintf('Test 1: Gaussian ...');
sample = exp(-linspace(-4,4,16).^2);
assert(all(abs(myfft(sample) - fft(sample))<1e-10));
fprintf('\tpassed\n');

fprintf('Test 1.1: Gaussian with a different Matlab dimension ...');
sample = exp(-linspace(-4,4,16).^2)';
assert(all(abs(myfft(sample) - fft(sample))<1e-10));
fprintf('\tpassed\n');

fprintf('Test 1.2: Gaussian complex ...');
sample = 1i*exp(-linspace(-4,4,128).^2);
assert(all(abs(myfft(sample) - fft(sample))<1e-10));
fprintf('\tpassed\n');

fprintf('Test 2: sawtooth ...');
sample = linspace(-1,1,128);
assert(all(abs(myfft(sample) - fft(sample))<1e-10));
fprintf('\tpassed\n');

fprintf('Test 3: sin and sin2 ...');
sample = sin(linspace(-pi,pi,128));
assert(all(abs(myfft(sample) - fft(sample))<1e-10));
sample = sin(2*linspace(-pi,pi,128));
assert(all(abs(myfft(sample) - fft(sample))<1e-10));
fprintf('\tpassed\n');

fprintf('All tests passed!\n');