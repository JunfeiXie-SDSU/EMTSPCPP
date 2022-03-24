function [globalMin,OptRoute,iter] = groupGA(pop,m,NoofRegions,object,Dmax,regions,...
    centralPs,InterRegionsPoints,numIter,popSize,Er,pm,TermCrit,term_method,time_max,toursearch)

% Find the Best pop in the Population
dists = vertcat(pop.totalDist);
[minDist,index] = min(dists);
globalMin = minDist;
OptRoute = pop(index).routes;



count = 0;
tic;
%% apply GA
for iter = 1:numIter
    numEr = floor((popSize - Er*popSize)/2)*2; 
    for k = 1:2:numEr
        %%%%%%%% Genetic Algorithm Operators
        %% parent selection
        indx_parents = GAselection(dists);

        %% crossover
        [child1,child2] = groupGAcrossover(indx_parents,pop,NoofRegions,m,object);
        %% mutation
        child1 = groupGAmutation(child1,m,pm);
        child2 = groupGAmutation(child2,m,pm);
        pop_new(k).gen =  child1.gen;
        pop_new(k+1).gen = child2.gen;
    end

    % evaluate the cost of each new solution
    for p = 1:numEr
        pop_new(p).totalDist = 0;
        pop_new(p).feasibility = ones(m,1);
        pop_new(p).routes = cell(m,1);
        pop_new(p).numarcs = zeros(m,1);
        pop_new(p).ltours = zeros(m,1);
        pop_new(p) = groupGAcost(pop_new(p),m,Dmax,object,regions,...
            centralPs,InterRegionsPoints,NoofRegions,iter,toursearch);
    end

    % Elitism
    [~, I] = sort(dists,'descend');
    pop(I(1:numEr)) = pop_new;
    
    %% evaluate the population
    dists = vertcat(pop.totalDist);
    [minDist,index] = min(dists);
    OptRoute = pop(index).routes;
    optFi = ~ismember(0,pop(index).feasibility);
    switch term_method 
        case 'count'
            if (globalMin - minDist) < 1e-1 && optFi == 1
                count = count + 1;
            else
                count = 0;
            end
            globalMin = minDist;

            if count >= TermCrit 
                break;
            end
        case 'time'
            time = toc; 
            if time > time_max
                break;
            end
            globalMin = minDist;
    end
end