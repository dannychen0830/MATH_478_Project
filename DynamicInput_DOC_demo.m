% Observe the effect of dynamic input on the manifestation of criticality.
% Here we have a cox process (gaussian white noise embedded in Poisson
% process) with random input drive and noise at each waiting interval

% Just press RUN

N = 40; 
W = zeros(N);
for i=1:N
    for j=1:N
        W(i,j) = 0.4*(rand - 0.5);
    end
end

dur = 500;
Vlast = zeros(N,1);
train = zeros([40,0]);
S = 0;

while true
    lambda = 1/(40 + normrnd(0,5));
    X = -1*log(1-rand)/lambda;
    S = S + X;
    drive = max([0.2+normrnd(0,0.3), 0.1]);
    noise = max([0.4 + normrnd(0,0.15), 0.2]);
    [spk NetParams V] = SimLIFNet(W,'simTime',X,'tstep',1e-2,...
      'offsetCurrents',drive*ones(length(W),1),...
           'noiseAmplitude',noise*ones(length(W),1),...
                'initialConditions',Vlast, ...
                'displayProgress',0,'plotResults',0);
    T = time_to_train(spk,X,1);
    if length(train(1,:)) == 0
        train = T;
    else
        train = [train(:,1:(length(train(1,:))-1)) train(:,length(train(1,:)))+T(:,1) T(:,2:length(T(1,:)))];
        train(train > 1) = 1;
    end
    Vlast = V(:,length(V));
    disp(['Epoch of length ', num2str(X)])
    if S > dur
        break 
    end
end

%%
figure(1)
subplot(2,3,[1 2 3])
imagesc(train)
title('spike train')
xlabel(''), ylabel('')
colormap(flipud(gray))
subplot(2,3,4)
[h,J] = InverseIsing2(train,N,30,0.4,0.4,10000,10000);
xlabel('number of iterations')
ylabel('discrepency between model and data')
title('Model Fitting')
subplot(2,3,6)
num = 5;
dis = zeros(1,num);
for i=1:num
    disp(['calculating distance: ',num2str(i)])
    dis(i) = heat_Capacity(h,J,N);
end
dim = [.85 0 .3 .2];
str = ['DOC: ',num2str(sum(dis)/num)];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
hold on
xline(1)
hold on 
%xline(sum(dis)/num + 1)
hold off
title('Temperature vs. Heat Capacity')
subplot(2,3,5)
J1 = J(:);
k = find(J1);
histogram(J1(k))
xlabel('pairwise parameter strength')
ylabel('frequency')
title('Distribution of Pairwise Parameter')
%%
figure(2)
samples = Metropolis_Ising(h,J,N,1,20000,20000);
s_mean = IsingMean(samples,N);
d_mean = DataMean(train);

[ha, Ja] = parseParam([(rand(1,N)-1.75)*(0.5) (rand(1,N^2-N)-0.5)*0.01]);
samples2 = Metropolis_Ising(ha,Ja,N,1,20000,20000);
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