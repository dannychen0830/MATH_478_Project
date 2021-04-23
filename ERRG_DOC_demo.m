N = 20;
eir = 0.75; % Excitatory/Inhibatory ratio
mag = 0.6; % Connection Strength
p = 0:0.1:0.7;

dist = zeros(1,length(p));
error = zeros(1,length(p));
for t = 1:length(p)
    num = 3;
    dis = zeros(1,num*5);
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
        [h,J] = InverseIsing(T,N,7,0.5,0.5,10000,10000);
        
        % calculate (average) distance to criticality
        %for i=1:num
        %    dis(i + num*(n-1)) = heat_Capacity(h,J,N);
        %end
        
        samples = Metropolis_Ising(h,J,N,1,30000,30000);
        energy = zeros(length(samples));
        for 
    end
    
    dist(t) = sum(dis)/(num*n);
    error(t) = std(dis);
    disp(['completed value ',num2str(p(t))])
end


errorbar(p,dist,error,'LineWidth',1.5)
xlabel('connection probability')
ylabel('distance to criticality')