% output the mean behavior of single & pairwise interactions 

% Input: samples - a vector of samples in binary representation

% Output: mean - a vector of mean and covariances of the sample spike train


function mean = IsingMean(samples,N)
ss = length(samples);
% convert to a spike train
train = num_to_spike(samples,N);
% convert to the {-1,1} formalism 
train(train==0) = -1;

% create a vector of means: <s_i>, <s_i s_j>
mean = zeros(1,N^2);
for i = 1:N
    mean(i) = sum(train(i,:))/ss;
end
index = N+1;
for i=1:N
    for j=i+1:N
        mean(index) = sum(train(i,:).*train(j,:))/ss - mean(i)*mean(j);
        index = index + 1;
    end
end