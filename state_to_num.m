% convert state vector to the corresponding integer 

% Input: x - a state vector in binary (-1,1) form

% Output: n - the same state in integer form

function n = state_to_num(x)
N = length(x);
x(x == -1) = 0;
n = 0;
for i = 1:N
    n = n + x(N-i+1)*2^(i-1);
end