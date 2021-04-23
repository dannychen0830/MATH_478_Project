% Processes the arrival times of independent Poisson Processes


function S = PoissonEpoch(trials,lambda,dur)

X = cell(trials,1); %interarrival
S = cell(trials,1); %arrival epoch
for i = 1:trials
    total = 0;
    j = 1;
    interarrival = zeros(1);
    while true
        xi = -1*log(1-rand)/lambda;
        total = total + xi;
        if total > dur
            break
        end
        interarrival(j) = xi;
        j = j + 1;
    end
    X(i) = {interarrival};
    S(i) = {cumsum(X{i})};
end