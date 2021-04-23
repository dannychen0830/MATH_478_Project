function mean = DataMean(data)

ss = length(data(1,:));
N = length(data(:,1));

data(data==0)=-1;

mean = zeros(1,N^2);
for i = 1:N
    mean(i) = sum(data(i,:))/ss;
end
index = N+1;
for i=1:N
    for j=i+1:N
        mean(index) = sum(data(i,:).*data(j,:))/ss - mean(i)*mean(j);
        index = index + 1;
    end
end