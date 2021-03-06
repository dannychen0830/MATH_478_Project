% Infer the parameters of Ising model (in the maximum entropy sense)
% given the spike train by minimizing the difference in the mean behavior

% Input: spike_train - a binary (0,1) representation of spiking, N - number
% of neurons, iter - the number of iteration desired, eta - learning rate
% (decrease at rate of O(1/t)), alpha - scaling of the inertia term, burnin
% - the length of the burn-in period for MCMC, ss - sample size of MCMC

% Output: h,J - the resulting parameter from the maximum entropy model

% Note: This is slightly different than InverseIsing in the sense that the
% parameters are updated imediately in each iteration. This should be the
% better algorithm though no significant difference was observed

function [h,J] = InverseIsing2(spike_train,N,iter,eta,alpha,burnin,ss)
g = [(rand(1,N)-1.75)*(0.5) (rand(1,N^2-N)-0.5)*0.01]; %vector of the coupling constants
inertia = zeros(1,N^2); %to carry the learning a bit further
data_mean = DataMean(spike_train); %the mean that the model is fitting

track_diff = zeros(1,iter);


for t=1:iter
    for i=1:N*(N+1)/2
        [h_t,J_t] = parseParam(g);
        samples = Metropolis_Ising(h_t,J_t,N,1,burnin,ss);
        iter_mean = IsingMean(samples,N);
        delta_g = -1*(eta/t)*(iter_mean(i) - data_mean(i)) + alpha*inertia(i);
        g(i) = g(i) + delta_g;
        inertia(i) = delta_g;
    end
    track_diff(t) = norm(iter_mean - data_mean); % track convergence
	%eta = eta/(t^(0.3));
    disp(['iteration ',num2str(t),' completed'])
end

[h,J] = parseParam(g);


plot(1:iter,track_diff,'linewidth',2)
xlabel('number of iterations')
ylabel('Data-to-Model Agreement')