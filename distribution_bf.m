% Analytically finds the Boltzmann distribution

% Inputs: N - number of neurons, Tk_b - temperature times boltzmann constant
% h - magnetic field, J - exchange coupling

% Outputs: P - the resulting boltzmann distribution

function P = distribution_bf(N,Tk_b,h,J)

if nargin < 3
    % randomly initialize 
    h = 2*rand(1,N) - 1;
    J = 2*randn(N,N) - 1;
    J = J - diag(diag(J)); % set diagonal as 0
    J = (J + transpose(J))/2; %make symmetric
end

% brute force calculate the Boltzmann distribution
P = zeros(2^N,1);
for i = 1:2^N
    P(i) = exp(hamiltonian(i-1,N,h,J)/Tk_b);
end
P = P/sum(P);

bar(0:2^N-1,P)