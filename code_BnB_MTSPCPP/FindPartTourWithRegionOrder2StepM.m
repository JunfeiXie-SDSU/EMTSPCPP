%%%% Given the visiting order of the regions, find the best intra-regional
%%%% path and entrance/exit locations using modified Two-Step
function [OptRoute, MinCost] = FindPartTourWithRegionOrder2StepM(regions, region_order, InterRegionsPoints, centralPs)
global UAS_initP

OptimalTT = region_order;
%% find the entrance and exit points at each region
% then find the inter-regional path
% initialize
OptRoute = UAS_initP;
distance = 0;
% find optimal tour
if OptimalTT(1) == 1
    for i = 2:length(OptimalTT)-1
        currRegion = OptimalTT(i);
        preRegion = OptimalTT(i-1);
        if preRegion == 1
            prePoint = UAS_initP;
        else
            prePoint = OptRoute(end,:);
        end
        nextRegion = OptimalTT(i+1);
        if nextRegion == 1
            nextPoint = UAS_initP;
        else
            nextPoint = centralPs(nextRegion-1,:);
        end
        numEEPairs = size(regions{currRegion-1}, 1)*4;
        currCost = zeros(numEEPairs, 1);
        for j = 1:numEEPairs
            tour = InterRegionsPoints{currRegion-1}{j};
            currCost(j) = norm(prePoint - tour(1,:)) + tour(end,1) + norm(tour(end-1,:) - nextPoint);
        end
        [nr,~] = find(currCost == min(currCost));
        OptRoute = [OptRoute; InterRegionsPoints{currRegion-1}{nr(1)}(1:end-1,:)];
        distance = distance + InterRegionsPoints{currRegion-1}{nr(1)}(end,1) + norm(prePoint - InterRegionsPoints{currRegion-1}{nr(1)}(1,:));
    end
    OptRoute = [OptRoute; UAS_initP];
    distance = distance + norm(OptRoute(end-1,:)-OptRoute(end,:));
    MinCost = distance;
else
    for i = 1:length(OptimalTT)-1
        currRegion = OptimalTT(i);
        preRegion = OptimalTT(end-1);
        if i == 1
            prePoint = centralPs(preRegion-1,:);
        else
            prePoint = OptRoute(end,:);
        end
       
        nextRegion = OptimalTT(i+1);
        if i ~= length(OptimalTT)-1
            nextPoint = centralPs(nextRegion-1,:);
        else
            nextPoint = OptRoute(1,:);
        end
        numEEPairs = size(regions{currRegion-1}, 1)*4;
        currCost = zeros(numEEPairs, 1);
        for j = 1:numEEPairs
            tour = InterRegionsPoints{currRegion-1}{j};
            currCost(j) = norm(prePoint - tour(1,:)) + tour(end,1) + norm(tour(end-1,:) - nextPoint);
        end
        [nr,~] = find(currCost == min(currCost));
        OptRoute = [OptRoute; InterRegionsPoints{currRegion-1}{nr(1)}(1:end-1,:)];
        if i == 1
            distance = distance + InterRegionsPoints{currRegion-1}{nr(1)}(end,1);
        else
            distance = distance + InterRegionsPoints{currRegion-1}{nr(1)}(end,1) + norm(prePoint - InterRegionsPoints{currRegion-1}{nr(1)}(1,:));
        end
    end
    OptRoute = [OptRoute; OptRoute(1,:)];
    distance = distance + norm(OptRoute(end-1,:)-OptRoute(end,:));
    MinCost = distance;
end
end