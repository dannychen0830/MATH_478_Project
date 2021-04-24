% Examine the how the signature of criticality changes with respect to
% different distributions of weights: Laplace, uniform, parabolic

% Just press RUN

% generate weight matrix (uncomment the corresponding line for the desired
% distribution)
N = 40;
W = zeros(N);
for i=1:N
    for j=1:N
        W(i,j) = Laplace(rand,0,0.1);%0.4*(rand - 0.5); %nthroot(rand-0.5,3)/4;
    end
end

% generate spike train
[spk NetParams V] = SimLIFNet(W,'simTime',500,'tstep',1e-2,...
      'offsetCurrents',0.1*ones(length(W),1),...
           'noiseAmplitude',0.6*ones(length(W),1));

T = time_to_train(spk,500,1);

% generate figures
figure(2)
subplot(1,4,1) % plot of weight distribution (histogram)
Wc = W(:);
histogram(Wc)
xlabel('weight between neurons')
ylabel('frequency')
title('Weight Distribution of Neural Network')

subplot(1,4,2) % plot of algorithm performance
[h,J] = InverseIsing2(T,N,100,0.1,0.1,8000,8000);
xlabel('number of iterations')
ylabel('discrepency between model and data')
title('Model Fitting')

subplot(1,4,4) % plot of heat capacity
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

subplot(1,4,3) % plot of distribution of pairwise parameter
J1 = J(:);
k = find(J1);
histogram(J1(k))
xlabel('pairwise parameter strength')
ylabel('frequency')
title('Distribution of Pairwise Parameter')
%% Take a closer look at the convergence
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


function s = Laplace(xi,mu,w)
if xi < 0.5
    s = mu + w*log(2*xi);
else
    s = mu - w*log(2*(1-xi));
end
end