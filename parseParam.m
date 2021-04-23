% sometimes the parameters are represented in vector form, this functoin
% writes it into a vector/matrix (h/J) so that it is compatible with the
% MCMC algorithm format

% Input: g - the parameters in one vector

% Output: h,J - the magnetic field and exchange couplings

function [h,J] = parseParam(g)
N = sqrt(length(g));
h = zeros(1,N);
J = zeros(N,N);
for i=1:N
    h(i) = g(i);
end
index = N+1;
for i=1:N
    for j=i+1:N
        J(i,j) = g(index);
        index = index + 1;
    end
end
