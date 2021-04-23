% takes a series of states and convert to spike train

% Input: num - vector of states in integer form, N - number of neurons

% Output: X - spike train of size N-by-(length of num) 

function X = num_to_spike(num,N)
s = length(num);
X = zeros(N,s);
% convert each column to state vector
for i=1:s
    X(:,i) = num_to_state(num(i),N);
end

% replace silent with 0
X(X==-1) = 0;
