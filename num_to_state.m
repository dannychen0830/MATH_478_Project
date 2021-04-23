% convert from an integer to the state vector

% Input: n - a state in integer form, N - number of neurons

% Output: x - the state in binary (-1, 1) form

function x = num_to_state(n,N)
x = zeros(N,1);
for i = N:-1:1
    r = rem(n,2);
    x(i) = r;
    n = (n-r)/2;
end

x(x==0) = -1;