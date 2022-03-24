%=========================================================================
% Branch-and-Bound method for EMTSP-CPP (Energy Constrainted Multiple TSP-CPP)
%
% Reference: J. Xie, J. Chen, "Multi-Regional Coverage Path Planning for Multiple Energy Constrained UAVs",
%           accepted by IEEE Transactions on Intelligent Transportation Systems, 2022.
%
% Date: 03/19/2022
%=========================================================================

%=========================================================================

%%%Functionality:
    % Generate the routes for multiple UAVs to cover multiple non-overlapping 
    % convex polygon regions 
    
%%%Input Dataset to be loaded:
    % 'nRegionsCase.mat' stores the map of the regions, where n is the
    % number of regions. The dataset includes 1) 'regions' that stores the
    % coordinates of the vertices of each region, and 2) 'NoofRegions' that
    % specifies the number of regions
%=========================================================================



clear all; clc;
load('ThreeRegionsCase.mat');
global UAS_range
global UAS_initP
UAS_range = [1.5,3];  % sensing range of the UAVs;
UAS_initP = [0,0]; %initial position of the UAVs. All UAVs depart from the same location;
m = 1;  % number of UAVs
Dmax = inf; % maximum distance each UAV can travel
searchmethod = 'DFS'; % Search strategy to use
% searchmethod = 'BFS'; % alternative search strategy for comparison
% purpose


InterRegionsPoints = FindAllBFPPaths3(NoofRegions, regions);

% find the central point for each region
centralPs = zeros(NoofRegions, 2);
for i = 1:NoofRegions
    rect = regions{i};
    centralPs(i,:) = mean(rect,1);
end

%% construct matrix D using the minimum distance between regions and for covering a region
% % find the minimum cost for covering a region
for i = 1:NoofRegions
    CandidatePaths = InterRegionsPoints{i};
    numPaths = length(CandidatePaths);
    CoveragePaths(i).npaths = numPaths;
    CoveragePaths(i).paths = cell(numPaths,1);
    CoveragePaths(i).pcosts = zeros(numPaths,1);
    CoveragePaths(i).startpoints = zeros(numPaths,2);
    CoveragePaths(i).endpoints = zeros(numPaths,2);
    for j = 1:length(CandidatePaths)
       CoveragePaths(i).paths{j} = CandidatePaths{j}(1:end-1,:);
       CoveragePaths(i).startpoints(j,:) = CandidatePaths{j}(1,:);
       CoveragePaths(i).endpoints(j,:) = CandidatePaths{j}(end-1,:);
       cost = CandidatePaths{j}(end,1);
       CoveragePaths(i).pcosts(j) = cost;
    end
end

% find the minimum distance between each region and the depot
D = zeros(NoofRegions+1, NoofRegions+1);
for i = 1:NoofRegions+1
    for j = 1:NoofRegions+1
       if i==j
           D(i,j) = inf;
       elseif i==1
           proRegion_indx = j-1;
           EndPoints = CoveragePaths(proRegion_indx).startpoints;
           dists = vecnorm(EndPoints - UAS_initP,2,2) + CoveragePaths(proRegion_indx).pcosts;
           D(i,j) = min(dists);
       elseif j == 1
           preRegion_indx = i-1;
           StartPoints = CoveragePaths(preRegion_indx).endpoints;
           dists = vecnorm(StartPoints - UAS_initP,2,2);
           D(i,j) = min(dists);
       else  % i~=j~=1
           preRegion_indx = i-1;
           proRegion_indx = j-1;
           StartPoints = CoveragePaths(preRegion_indx).endpoints;
           EndPoints = CoveragePaths(proRegion_indx).startpoints;
           [pre, pro] = ndgrid(1:CoveragePaths(preRegion_indx).npaths, 1:CoveragePaths(proRegion_indx).npaths);
           comb = [pre(:),pro(:)];
           dis_comb = vecnorm(StartPoints(comb(:,1),:) - EndPoints(comb(:,2),:),2,2) + ...
               CoveragePaths(proRegion_indx).pcosts(comb(:,2));
           D(i,j) = min(dis_comb);
       end
    end
end

Di = D;
% create new distance matrix by adding m-1 copy of the depot
for i = 1:m-1
    D = [Di, Di(:,1)]; %add column
    D = [D; Di(1,:),inf]; %add row
    Di = D;
end

% root node is the depot
ii = 1; % index of the tree
[cost, reducedD] = calculateCost(D);
STree(ii).vertex = 1; % stores current region index
STree(ii).level = 1; % number of regions visited so far
STree(ii).reducedMatrix = reducedD;
STree(ii).LB = cost; % store the lower bound
STree(ii).tour = 1; % tour constructed so far
STree(ii).route = 1; % index of path in each region
UB = m*Dmax;
iter = 1;
while ~isempty(STree) 
    iter = iter + 1;
    switch searchmethod
        case 'BFS'
            %% find a live node bi with the least estimated cost
            L_LB = vertcat(STree.LB);
            Cmin = min(L_LB);
            indx = find(L_LB == Cmin);
            i = randi(length(indx),1,1);
            bi = indx(i);
            currRegion = STree(bi).vertex;% index of the current region + m
            %% if all cities are visited
            if STree(bi).level == NoofRegions+m
                OptRoute_indx = STree(bi).route;
                OptTour = STree(bi).tour;
                OptLB = STree(bi).LB;
                break;
            end
        case 'DFS'
            %% find a live node bi with the maximum number of levels
            L_level = vertcat(STree.level);
            L_LB = vertcat(STree.LB);
            indx = find(L_level == max(L_level));
            if length(indx)>1
                % pick the one with the lowest LB
                L_LB_indx = L_LB(indx);
                indxx = find(L_LB_indx == min(L_LB(indx)));
                bi = indx(indxx(1));
            else
                bi = indx;
            end
            currRegion = STree(bi).vertex;% index of the current region + m
            if STree(bi).level == NoofRegions+m
                if STree(bi).LB <= UB
                    OptRoute_indx = STree(bi).route;
                    OptTour = STree(bi).tour;
                    OptLB = STree(bi).LB;
                end
                STree(bi) = [];
                ii = length(STree);
                continue;
            end
                
    end
    currRoute = STree(bi).route;
    update = 0;
    regions_left = 1:NoofRegions+m;
    regions_left(STree(bi).tour) = [];
    noofregions_left = length(regions_left);
    costs = zeros(noofregions_left,1);
    %% expand node bi
    for j = 1:noofregions_left
        nextRegion = regions_left(j);
        if (STree(bi).reducedMatrix(currRegion,nextRegion)) ~= inf
            [costs(j),~] = costIncrease(STree(bi).reducedMatrix, currRegion, nextRegion);
        else
            costs(j) = -1;
        end
    end
    costInx = find(costs == max(costs));
    if ~isempty(costInx) && max(costs) > -1
        if length(costInx) > 1
            i = randi(length(costInx),1,1);
            j = regions_left(costInx(i));
        else
            j = regions_left(costInx);
        end
        ii = ii + 1;
        STree(ii) = STree(bi);
        [~, STree(ii).reducedMatrix] = costIncrease(STree(ii).reducedMatrix, currRegion, j);
        STree(ii).LB = STree(ii).LB + max(costs);
        if STree(ii).LB > UB
            STree(ii) = [];
            ii = ii - 1;
        end
            % create a child node 
            ii = ii + 1;
            STree(ii) = newNode(STree(bi),j);
            [cost, STree(ii).reducedMatrix] = calculateCost(STree(ii).reducedMatrix);
            %% update the LB
            if currRegion == 1
                endpoint = UAS_initP;
                % pick the path with the shortest distance
                startpoints = CoveragePaths(j-1).startpoints;
                dists = vecnorm(endpoint - startpoints,2,2)+CoveragePaths(j-1).pcosts;
                [dmin, dmin_indx] = min(dists);
                % calculate the lower bound of the path starting at node j
                STree(ii).LB = STree(bi).LB + cost+...
                    dmin - D(currRegion,j) + STree(bi).reducedMatrix(currRegion,j);
                STree(ii).route = [STree(ii).route,dmin_indx];
            else
                preRegion = STree(bi).tour(end-1);
                if preRegion == 1 || preRegion > NoofRegions + 1
                    prePoint = UAS_initP;
                else
                    prePoint_pathindx = STree(bi).route(end-1);
                    prePoint = CoveragePaths(preRegion-1).endpoints(prePoint_pathindx,:);
                end
                if j <= NoofRegions + 1
                    endPoint = centralPs(j-1,:);
                else
                    endPoint = UAS_initP;
                end
                if currRegion <= NoofRegions + 1
                    sPoints = CoveragePaths(currRegion-1).startpoints;
                    ePoints = CoveragePaths(currRegion-1).endpoints;
                    dists = vecnorm(prePoint - sPoints,2,2) + CoveragePaths(currRegion-1).pcosts...
                        + vecnorm(endPoint - ePoints,2,2);
                    [~, dmin_indx] = min(dists);
                    dmin = vecnorm(prePoint - sPoints(dmin_indx,:),2,2) + CoveragePaths(currRegion-1).pcosts(dmin_indx);
                    currPath_indx = STree(bi).route(end);
                    currStartpoint = CoveragePaths(currRegion-1).startpoints(currPath_indx,:);
                    deltaD = dmin-(vecnorm(prePoint-currStartpoint,2,2)+CoveragePaths(currRegion-1).pcosts(currPath_indx));
                    currLB = STree(bi).LB + deltaD;
                    endpoint = CoveragePaths(currRegion-1).endpoints(dmin_indx,:);
                else
                    currLB = STree(bi).LB;
                    endpoint = UAS_initP;
                end
                %% pick the path with the shortest distance
                if j <= NoofRegions + 1
                    startpoints = CoveragePaths(j-1).startpoints;
                    if STree(bi).level < NoofRegions+m-1
                        dists_new = vecnorm(endpoint - startpoints,2,2)+CoveragePaths(j-1).pcosts;
                        [dmin_new, dmew_indx] = min(dists_new);
                        STree(ii).LB = currLB + cost+...
                        dmin_new - D(currRegion,j) + STree(bi).reducedMatrix(currRegion,j);
                        % check feasibility
                        [~,tindx] = find(STree(ii).tour == 1 | STree(ii).tour > NoofRegions+1);
                        startr = tindx(end)+1;
                        subtour = [1, STree(ii).tour(startr:end),1];
                        [~, c] = FindPartTourWithRegionOrder2StepM(regions,...
                            subtour, InterRegionsPoints, centralPs);
                        if c > Dmax
                            STree(ii).LB = inf;
                        end
                    else %% last level, which leads to a complete tour
                        dists_new = vecnorm(endpoint - startpoints,2,2)+CoveragePaths(j-1).pcosts...
                            + vecnorm(UAS_initP - CoveragePaths(j-1).endpoints,2,2);
                        [dmin_new, dmew_indx] = min(dists_new);
                        node_temp = newNode(STree(ii),1);
                        STree(ii).LB = currLB + cost + ...
                        dmin_new - D(currRegion,j)-D(j,1)+ STree(bi).reducedMatrix(currRegion,j);
                        % check feasibility
                        [~,tindx] = find(STree(ii).tour == 1 | STree(ii).tour > NoofRegions+1);
                        startr = tindx(end)+1;
                        subtour = [1, STree(ii).tour(startr:end),1];
                        [~, c] = FindPartTourWithRegionOrder2StepM(regions,...
                            subtour, InterRegionsPoints, centralPs);
                        if c > Dmax
                            STree(ii).LB = inf;
                        end
                        if STree(ii).LB < UB
                            UB = STree(ii).LB;
                            update = 1;
                        end
                    end
                else
                    if STree(bi).level < NoofRegions + m - 1
                        dmin_new = vecnorm(endpoint - UAS_initP,2,2);
                        STree(ii).LB = currLB + cost + dmin_new - D(currRegion,j)...
                            +STree(bi).reducedMatrix(currRegion,j);
                        dmew_indx = 1;
                        %% this indicates the end of one tour and the start of another tour
                        [~,tindx] = find(STree(ii).tour == 1 | STree(ii).tour > NoofRegions+1);
                        % check feasibility
                        startr = tindx(end-1)+1;
                        endr = tindx(end)-1;
                        subtour = [1, STree(ii).tour(startr:endr),1];
                        [~, c] = FindPartTourWithRegionOrder2StepM(regions,...
                            subtour, InterRegionsPoints, centralPs);
                        if c > Dmax
                            STree(ii).LB = inf;
                        end
                    else
                        STree(ii).LB = inf;
                        dmew_index = 1;
                    end
                end
                STree(ii).route = [STree(ii).route(1:end-1),dmin_indx,dmew_indx];
            end

            if STree(ii).LB > UB
                STree(ii) = [];
                ii = ii - 1;
            end

    end
    if update == 1 
       lindx = find(L_LB > UB);
       STree(unique([lindx;bi])) = [];
    else
        STree(bi) = [];
    end
    ii = length(STree);
    
end


if UB == m*Dmax
    ind = 2;
else
    OptRoute = UAS_initP;
    for i = 2:NoofRegions+m
        currR = OptTour(i);
        if currR <= NoofRegions + 1
            path = CoveragePaths(currR-1).paths{OptRoute_indx(i)};
        else
            path = [UAS_initP;UAS_initP];
        end
        OptRoute = [OptRoute;path];
    end
    OptRoute = [OptRoute;UAS_initP];
    ind = 1;
end


colors =  ['b','k','y','c','g','r','m'];
figure;
hold on;
for i = 1:NoofRegions
    rectangle = regions{i};
    fill(rectangle(:,1),rectangle(:,2),[0.9,0.9,0.9])
end

ltours = zeros(m,1);
l = 1;
for i = 1:length(OptRoute(:,1))-1
    if ~isequal(OptRoute(i,:), OptRoute(i+1,:))
        line(OptRoute(i:i+1,1), OptRoute(i:i+1,2), 'Color', colors(l), 'Marker','.', 'LineWidth', 1.5);
        ltours(l) = ltours(l) + norm(OptRoute(i,:)-OptRoute(i+1,:));
    else
        l = l+1;
    end
end
plot(OptRoute(1,1), OptRoute(1,2), '>', 'MarkerSize',10, 'MarkerFaceColor', 'r');
hold off;

totalcost = sum(ltours)
ltours


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function node = newNode(parent,j)
node.tour = [parent.tour,j];
node.level = parent.level+1;
node.LB = parent.LB;
node.route = parent.route;
node.reducedMatrix = parent.reducedMatrix;
N = size(node.reducedMatrix,1);
i = parent.vertex;  % node j, travel from i
% change all entries of row i and column j to infinity
for k = 1:N
    node.reducedMatrix(i,k) = inf;
    node.reducedMatrix(k,j) = inf;
end
% node.reducedMatrix(j,1) = inf;
node.reducedMatrix(j,i) = inf;
node.vertex = j;
end

function [reducedMatrix,rmin] = rowReduction(reducedMatrix)
N = size(reducedMatrix,1);
% find minimum value of each row
rmin = min(reducedMatrix,[],2);
% reduce the minimum value from each element in each row
for i = 1:N
    for j = 1:N
        if reducedMatrix(i,j)~=inf && rmin(i)~=inf
            reducedMatrix(i,j) = reducedMatrix(i,j) - rmin(i); 
        end
    end
end
end


function [reducedMatrix,cmin] = columnReduction(reducedMatrix)
N = size(reducedMatrix,1);
% find minimum value of each row
cmin = min(reducedMatrix,[],1);
% reduce the minimum value from each element in each row
for i = 1:N
    for j = 1:N
        if reducedMatrix(i,j)~=inf && cmin(j)~=inf
            reducedMatrix(i,j) = reducedMatrix(i,j) - cmin(j); 
        end
    end
end
end

% function to get the lower bound
function [cost, reducedMatrix] = calculateCost(reducedMatrix)

[reducedMatrix,rmin] = rowReduction(reducedMatrix);
[reducedMatrix,cmin] = columnReduction(reducedMatrix);
cost = sum(rmin(rmin~=inf)) + sum(cmin(cmin~=inf));
end

function [cost, reducedMatrix] = costIncrease(reducedMatrix, i, j)
reducedMatrix(i,j) = inf;
[cost, reducedMatrix] = calculateCost(reducedMatrix);
end


