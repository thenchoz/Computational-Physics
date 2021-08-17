function result = mydft(input_array)
% mydft : computes discrete Fourier transform.
%
% Arguments:
%       input array (1D array complex): data to transform;
%
% Returns : 1D array complex,
%       transformed data of the same shape as an input array.


N = length(input_array);
result = input_array; % having right size

for m = 1:N
  result(m) = 0;
  for n = 1:N
    % m, n -1 cause of matlab index
    result(m) = result(m) + input_array(n) * exp(-2*pi*1i * (m-1)*(n-1)/N);
  end
end
end
