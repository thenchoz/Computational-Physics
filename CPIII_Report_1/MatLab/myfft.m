function result=myfft(input_array)
% myfft : computes discrete Fourier transform with radix-2.
%
% Arguments:
%       input array (1D array complex, size 2^N, N>0): data to transform.
%
% Returns : 1D array complex,
%       transformed data of the same shape as an input array.


N = length(input_array);

if(N >= 2)
    % splitting array in two
    even = input_array(1:2:N);
    odd  = input_array(2:2:N);
    
    % applying recursion
    even = myfft(even);
    odd  = myfft(odd);
    
    % fourier transform
    for j = 1:N/2
        input_array(j)       = even(j) + W(j-1, N) * odd(j);
        input_array(j + N/2) = even(j) - W(j-1, N) * odd(j);
    end
end
   result = input_array;
end

function result = W(k, N)
% W : exp as in lecture notes

result = exp(- 2 * pi * 1i * k / N);
end
