% Uses the Metropolis-Hastings algorithm to sample Ising model

% Input: h - magnetic field, J - exchange couplings, N - number of neurons,
% TK-b - temperature times boltzmann constant, nstep - the length of the
% burn-in period, nsample - the number of samples desired

function samples = Metropolis_Ising(h,J,N,Tk_b,nstep,nsample)
% pick an inital state
x_0 = zeros(N,1);
for i = 1:N
    if rand <= 0.5, x_0(i) = 1;
    else, x_0(i) = -1;
    end
end

samples = zeros(1,nsample);

% MIX THE CHAIN! (for nstep) + record outcome
for i = 1:nstep + nsample
    % find proposal state
    x_1 = x_0;
    r = randi(N);
    x_1(r) = -1*x_1(r);
    %for j = 1:N
     %   if rand <= 0.5, x_1(j) = 1;
      %  else, x_1(j) = -1;
       % end
    %end
    
    % find transition probability
    dE = hamiltonian(x_1,N,h,J) - hamiltonian(x_0,N,h,J);
    transition_prob = min(1,exp(dE/Tk_b));
    
    if rand <= transition_prob, x_0 = x_1; end
    %if state_to_num(x_0) == 0 && 0 ~= state_to_num(x_1), count = count + 1; end
    % when it's recording time, record the sample
    if i > nstep, samples(i-nstep) = state_to_num(x_0); end
end 
