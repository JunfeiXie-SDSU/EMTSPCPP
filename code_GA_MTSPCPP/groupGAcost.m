function pop = groupGAcost(pop,m,Dmax,object,regions,centralPs,InterRegionsPoints,...
    NoofRegions,iter,toursearch)

%% use for objective 1, dynamic penalty
c_dynamic = 0.5;
alpha_dynamic = 1;
beta_dynamic = 2;


for i = 1:m
    indx = pop.gen{i};
    [region_order, d, route] = findSubtour2(indx, regions, ...
    InterRegionsPoints,centralPs,NoofRegions,toursearch);
    pop.ltours(i) = d;
    pop.gen{i} = region_order(2:end-1) - ones(1,length(indx));
    pop.routes{i} = route;
    pop.numarcs(i) = length(indx);
    if object == 1
    % objective function 1
        pop.totalDist = pop.totalDist + d;
        if Dmax < inf && d > Dmax
            pop.totalDist = pop.totalDist + (c_dynamic*iter)^alpha_dynamic*(max(0, (d - Dmax)))^beta_dynamic;
        end
    else
    % objective function 2
        pop.totalDist = max(pop.totalDist,d);
    end
    if d > Dmax
        pop.feasibility(i) = 0;
    end
end
    
end