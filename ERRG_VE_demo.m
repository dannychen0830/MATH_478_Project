% Plots the variance of energy with respect to network density (governed by
% the connection probability in a ER random graph)

% Just press RUN to see the plot of variance vs. p


N = 20;
eir = 0.75; % Excitatory/Inhibatory ratio
mag = 0.6; % Connection Strength
p = 0:0.1:0.7;

VE = zeros(1,length(p));
error = zeros(1,length(p));
for t = 1:length(p)
    ss = 15000;
    energy = zeros(5,ss);
    for n=1:5
        % create ER random graph
        W = zeros(N);
        for i = 1:N
            for j = i+1:N
                if rand < p(t)
                    if rand < eir
                        if rand < 0.5, W(i,j) = mag; else, W(j,i) = mag; end
                    else
                        if rand < 0.5, W(i,j) = -1*mag;else, W(j,i) = -1*mag; end
                    end
                end
            end
        end
    
        % simulate spikes
        [spk NetParams V] = SimLIFNet(W,'simTime',100,'tstep',1e-2,...
          'offsetCurrents',0.1*ones(length(W),1),...
          'noiseAmplitude',0.6*ones(length(W),1),...
          'displayProgress',0,'plotResults',0);
        T = time_to_train(spk,100,1);
        
        % fit Ising model to spikes
        [h,J] = InverseIsing(T,N,50,0.1,0.1,5000,5000);
        
        samples = Metropolis_Ising(h,J,N,1,10000,ss);
        for k=1:ss
            energy(n,:) = hamiltonian(samples(k),N,h,J);
        end
    end
    
    VE(t) = var(energy(:));
    error(t) = std(var(energy));
    disp(['completed value ',num2str(p(t))])
end

errorbar(p,VE,error,'LineWidth',1.5)
xlabel('connection probability')
ylabel('variance of energy')