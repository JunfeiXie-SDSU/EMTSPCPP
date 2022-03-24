function [indx_parents] = GAselection(totalDist)
popSize = length(totalDist);

[~,ranks] = sort(totalDist);
fitness = popSize:-1:1;
total_fitness = sum(fitness);
indx_parents = zeros(2,1);
p = 1;
while p <= 2
    num = total_fitness * rand(1,1);
    P = 0;
    i = 0;
    while P < num
        i = i+1;
        P = fitness(i) + P;
    end
    indx_parents(p) = ranks(i);
    if length(unique(indx_parents(1:p))) == p 
        p = p+1;
    end
end
