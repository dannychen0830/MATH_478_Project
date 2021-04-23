% Comparison between MCMC and actual distribution with randomized
% parameters. The whole distribution is plotted and we can see that MCMC
% captures the important "features" of the distribution.

function MCMC_demo

N=8; Tk_b=4;

% randomly initialize 
h = 2*rand(1,N);
J = zeros(N,N);
for i=1:N
    for j=i+1:N
        J(i,j) = 2*rand - 1;
    end
end

subplot(1,2,1)
distribution_bf(N,Tk_b,h,J)
title('Actual distribution')
xlabel('state')
xlim([0 2^N])


sample = Metropolis_Ising(h,J,N,Tk_b,10000,10000);
subplot(1,2,2)
histogram(sample,2^N)
xlim([0 2^N])
title('Metropolis-Hastings')
xlabel('state')