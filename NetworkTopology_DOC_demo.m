% Examine the how the signature of criticality changes with respect to
% different network topology: Erdos-Renyi Random Graph, Small World Network

% Go to the desired section for creating different graphs, then go to the
% simulation section


%% Erdos-Renyi Random Graph
N = 40;
p = 0.15;
eir = 0.7;
mag = 0.45;

A = zeros(N);
W = zeros(N);

for i=1:N
    for j=i+1:N
        if rand < p
            A(i,j) = 1; A(j,i) = 1;
            if rand < eir
                if rand < 0.5, W(i,j) = mag; else, W(j,i) = mag; end
            else 
                if rand < 0.5, W(i,j) = -mag; else, W(j,i) = mag; end
            end
        end
    end
end

G = graph(A);
subplot(1,2,1)
plot(G,'NodeLabel',{})
title('Resulting Network')
subplot(1,2,2)
histogram(degree(G))
title('Degree distribution') % note that for small Ns, this might not look binomial

%% Small World Network

N = 40;
p1 = 0.65;
p2 = 0.02;
adj = 3;
eir = 0.70;
mag = 0.45;

A = zeros(N);
W = zeros(N);

% connect adjacent edges
for i = 0:N-1
    for j = 1:adj
        if rand < p1 
            k = rem(i + j,40) + 1; 
            A(i+1,k) = 1; A(k,i+1) = 1;
            if rand < eir
                if rand < 0.5, W(i+1,k) = mag; else, W(k,i+1) = mag; end
            else 
                if rand < 0.5, W(i+1,k) = -mag; else, W(k,i+1) = mag; end
            end
        end
    end
end

% connect clusters
for i = 1:N
    for j=i+1:N
        if rand < p2
            A(i,j) = 1; A(j,i) = 1;
            if rand < eir
                if rand < 0.5, W(i,j) = mag; else, W(j,i) = mag; end
            else 
                if rand < 0.5, W(i,j) = -mag; else, W(j,i) = mag; end
            end
        end
    end
end

G = graph(A);
subplot(1,2,1)
plot(G,'NodeLabel',{})
title('Resulting Network')
subplot(1,2,2)
histogram(degree(G))
title('Degree distribution')

%% Simulation 
[spk NetParams V] = SimLIFNet(W,'simTime',500,'tstep',1e-2,...
      'offsetCurrents',0.1*ones(length(W),1),...
           'noiseAmplitude',0.6*ones(length(W),1));
%%
T = time_to_train(spk,500,1);
figure(2)
subplot(1,4,1)
plot(G,'NodeLabel',{})
title('Network Topolgy')
subplot(1,4,2)
[h,J] = InverseIsing2(T,N,100,0.1,0.1,20000,20000);
xlabel('number of iterations')
ylabel('discrepency between model and data')
title('Model Fitting')
subplot(1,4,4)
num = 5;
dis = zeros(1,num);
for i=1:num
    disp(['calculating distance: ',num2str(i)])
    dis(i) = heat_Capacity(h,J,N);
end
dim = [.83 0 .3 .3];
str = ['DOC: ', num2str(sum(dis)/num)];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
hold on
xline(1)
hold on 
xline(sum(dis)/num + 1)
hold off
title('Temperature vs. Heat Capacity')
subplot(1,4,3)
J1 = J(:);
k = find(J1);
histogram(J1(k))
xlabel('pairwise parameter strength')
ylabel('frequency')
title('Distribution of Pairwise Parameter')
%% Examine Convergence
figure(3)
samples = Metropolis_Ising(h,J,N,1,10000,10000);
s_mean = IsingMean(samples,N);
d_mean = DataMean(T);

[ha, Ja] = parseParam((rand(1,N^2)-0.5)*0.01);
samples2 = Metropolis_Ising(ha,Ja,N,1,10000,10000);
s_mean2 = IsingMean(samples2,N);

plot(d_mean(N+1:N*(N+1)/2),s_mean2(N+1:N*(N+1)/2),'b.','MarkerSize',8)
hold on
plot(d_mean(N+1:N*(N+1)/2),s_mean(N+1:N*(N+1)/2),'r.','MarkerSize',8)
hold on
plot(-0.1:0.01:0.1,-0.1:0.01:0.1)
hold off
xlabel('data correlation')
ylabel('model correlation')
legend('initial','fitted')
