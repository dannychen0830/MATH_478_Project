% This function converts a array of cells containing spike times into spike
% trains

% Input: spk - an array of cells containing the spike time for each cell,
% dur - the duration of the recording,  tbin - the length of each time bin

function train = time_to_train(spk,dur,tbin)

nbin = ceil(dur/tbin);
train = zeros(length(spk),nbin);
for i=1:length(spk)
    for j=1:length(spk{i})
        ind = ceil(spk{i}(j)/tbin);
        train(i,ind) = train(i,ind) + 1;
    end
end

if ~ isempty(find(train>=2,1))
    disp("Multiple Spikes! Check or Reduce Bin Size!")
    train(train>=2) = 1;
end