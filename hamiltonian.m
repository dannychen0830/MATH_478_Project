% calculate the hamiltonian of pairwise Ising model

% Inputs: n - the state vector (in either binary or integer
% representation), N - number of neurons, h - magnetic field, J - exchange
% couplings

% Output: H - the hamiltonian, or energy, of the state

function H = hamiltonian(n, N, h, J)

% if input is in integer representation, making it into state
if length(n) < 2
    x = num_to_state(n,N); 
% otherwise, keep it in state form
else
    x = n;
end

% calculate the hamiltonian
H = h*x + transpose(x)*J*x;