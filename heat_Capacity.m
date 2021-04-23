% Plots the "Signature of Criticality" by varying the temperature and
% plotting the corresponding heat capacity (variance in energy)

% Input: h,J - magnetic field and exchange couplings, N - number of neurons

% Output: Plot of heat capacity and a polynomial fit to the data points for
% visualization purposes

% Note that we just approximated the function with a 10-degree polynomial,
% there is no rigorous justification for it (Weierstrass Approximation
% theorem?), it is just for visualization...

function [distance] = heat_Capacity(h,J,N)

dt = 0.1;
t_range = 0.1:dt:4;
C = zeros(1,length(t_range));

for t = 1:length(t_range)
    samples = Metropolis_Ising(h,J,N,t_range(t),20000,20000);
    E1 = zeros(1,length(samples));
    E2 = zeros(1,length(samples));
    for i=1:length(samples)
        E1(i) = hamiltonian(samples(i),N,h,J);
        E2(i) = hamiltonian(samples(i),N,h,J)^2;
    end
    C(t) = (sum(E2)/length(E2) - (sum(E1)/length(E1))^2)/t_range(t);
end

f = polyfit(t_range,C,10);
curve = polyval(f,t_range);
plot(t_range,C,'.')
hold on
plot(t_range,curve,'linewidth',1.5)
hold off

xlabel('temperature')
ylabel('heat capacity')

[val, maxind] = max(curve);
m = t_range(maxind);

[cs, ind] = sort(C);
last = ind((length(t_range)-4:length(t_range)));
last = t_range(last);
distance = (0.6*m + 0.4*sum(last)/5) - 1;