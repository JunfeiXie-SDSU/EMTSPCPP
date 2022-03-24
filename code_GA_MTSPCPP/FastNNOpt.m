function [OptRoute, distance, order_indx] = FastNNOpt(NoofRegions, regions, InterRegionsPoints,centralPs)
global UAS_initP

grids = [UAS_initP; centralPs];

[OptimalT,~] = tsp_nn(grids);
[~, c] = find(OptimalT == 1);
OptimalTT = [OptimalT(c:end), OptimalT(1:c-1),1];
% calculate distance table
cities = [UAS_initP; centralPs];
Dtable = zeros(NoofRegions+1, NoofRegions+1);
for i = 1:1+NoofRegions
    dists = vecnorm(kron(cities(i,:), ones(NoofRegions+1,1)) - cities,2,2);
    Dtable(i,:) = dists;
end

%% find the TSP-CPP tour
region_order = OptimalTT;
[OptRoute, dmin] = FindTourWithRegionOrder2StepM(regions, region_order, InterRegionsPoints,centralPs);
n = NoofRegions + 1;
p = OptimalTT(1:end-1);
cost = Inf;
pmin = p;
while cost - dmin > 1e-2
    cost = dmin;
    i = 0;
    b = p(n);
    % Loop over all edge pairs (ab,cd)
    while i < n-2
        a = b;
        i = i+1;
        b = p(i);
        j = i+1;
        d = p(j);
        while j < n
            c = d;
            j = j+1;
            d = p(j);
            % if (dist(ac) < dist(ab))
            dist1 = Dtable(a,c);
            dist2 = Dtable(a,b);
            if dist1 < dist2
                p_temp = p;
                % alternative path
                p_temp(i:j-1) = p_temp(j-1:-1:i);
                if i == 1
                    p_temp = [p_temp(j-1:end), p_temp(1:j-2)];
                end
                region_order = [p_temp,1];
                [Tour_temp, d_temp] = FindTourWithRegionOrder2StepM(regions, region_order, InterRegionsPoints,centralPs);
                % Keep best exchange
                if d_temp < dmin
                    dmin = d_temp;
                    OptRoute = Tour_temp;
                    pmin = p_temp;
                end
            end
        end
    end
    p = pmin;
end

distance = 0;
for i = 1:length(OptRoute(:,1))-1
    if i >= 2
        distance = distance + norm(OptRoute(i,:)-OptRoute(i+1,:));
    else
        distance = distance + norm(OptRoute(i,:)-OptRoute(i+1,:));
    end
end

order_indx = pmin(2:end) - ones(1, length(pmin)-1);
