function pop = InitGroup(m, NoofRegions, popSize, Dmax, regions,InterRegionsPoints,centralPs,object,toursearch)
global UAS_initP
method = 2;
%% construct matrix D using the minimum distance between regions and for covering a region
% % find the minimum cost for covering a region
for i = 1:NoofRegions
    CandidatePaths = InterRegionsPoints{i};
    numPaths = length(CandidatePaths);
    CoveragePaths(i).npaths = numPaths;
    CoveragePaths(i).pcosts = zeros(numPaths,1);
    CoveragePaths(i).startpoints = zeros(numPaths,2);
    CoveragePaths(i).endpoints = zeros(numPaths,2);
    for j = 1:length(CandidatePaths)
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

[theta,~] = cart2pol(centralPs(:,1)-UAS_initP(1), centralPs(:,2)-UAS_initP(2));
X = [cos(theta), sin(theta)];

pop = struct;

dmin = min(D(:,2:end),[],'all');

if method == 1 % suitable for object 1
    k = 1;
    while k <= popSize 
        unvisited = 1:NoofRegions;
        gene = cell(m,1);
        dist = zeros(m,1);
        while ~isempty(unvisited)
            % randomly pick a region
            ri = randi(length(unvisited));
            r = unvisited(ri);
            if object == 1
                % pick the salesmen that leads to the smallest increase of the
                % total distance
                delta_d = zeros(m,1);
                for i = 1:m
                    indx = [gene{i},r];
                    [~, d, ~] = findSubtour2(indx, regions, ...
                    InterRegionsPoints,centralPs,NoofRegions,'nn');
                    if d > Dmax
                        delta_d(i) = inf;
                    else
                        delta_d(i) = d - dist(i);
                    end
                end
                [~,minx] = min(delta_d);
                gene{minx} = [gene{minx},r];
                dist(minx) = dist(minx) + delta_d(minx);
            else
                % pick the salesmen that leads to the minimum longest path
                delta_d = zeros(m,1);
                dist_d = zeros(m,1);
                for i = 1:m
                    indx = [gene{i},r];
                    [~, d, ~] = findSubtour2(indx, regions, ...
                    InterRegionsPoints,centralPs,NoofRegions,'nn');
                    if d > Dmax
                        delta_d(i) = inf;
                    else
                        delta_d(i) = max([dist;d]);
                        dist_d(i) = d;
                    end
                end
                [~,minx] = min(delta_d);
                gene{minx} = [gene{minx},r];
                dist(minx) = dist_d(minx);
            end
            unvisited = unvisited(~ismember(unvisited,r));
        end
        pop(k).gen = gene;
        pop(k).totalDist = 0;
        pop(k).feasibility = ones(m,1);
        pop(k).routes = cell(m,1);
        pop(k).numarcs = zeros(m,1);
        pop(k).ltours = zeros(m,1);
        pop(k) = groupGAcost(pop(k),m,Dmax,object,regions,centralPs,InterRegionsPoints,NoofRegions,1,toursearch);
        if ~ismember(0,pop(k).numarcs)
            k = k+1;
        end
    end 
else
    k = 1;
    while k <= popSize % randomly pick k regions
        unvisited = 1:NoofRegions;
        pop(k).gen = cell(m,1);
        pop(k).ltours = zeros(m,1);
        centrals = randperm(NoofRegions, m);
        unvisited = unvisited(~ismember(unvisited,centrals));
        unassigned = 1:m;
        for i = 1:m
            pop(k).gen{i} = [pop(k).gen{i},centrals(i)];
            pop(k).ltours(i) = pop(k).ltours(i) + D(1,centrals(i)+1);
            if pop(k).ltours(i)+D(centrals(i)+1,1) + dmin > Dmax
                unassigned = unassigned(~ismember(unassigned,i));
            end
        end
        
        if object == 1
            while ~isempty(unvisited)
                if ~isempty(unassigned)
                    dist_m = zeros(length(unassigned),2);
                    for l = 1:length(unassigned)
                        i = unassigned(l);
                        dists = vecnorm(X(unvisited,:)-X(centrals(i),:),2,2);
                        [md,indx] = min(dists);
                        r = unvisited(indx);
                        dist_m(l,:) = [md,r];
                    end
                    [~,minx] = min(dist_m(:,1));
                    s = unassigned(minx);
                    r = dist_m(minx,2);
                    pre = pop(k).gen{s}(end);
                    if pop(k).ltours(s)+D(pre+1,r+1) + dmin > Dmax
                        unassigned = unassigned(~ismember(unassigned,s));
                    else               
                        pop(k).gen{s} = [pop(k).gen{s},r];
                        pop(k).ltours(s) = pop(k).ltours(s) + D(pre+1,r+1);
                        unvisited = unvisited(~ismember(unvisited,r));
                        if isempty(unvisited) 
                            break;
                        end
                    end
                else
                    len_unvisited = length(unvisited);
                    for i = 1:len_unvisited
                        r = unvisited(i);
                        indx = randi(m);
                        pop(k).gen{indx} = [pop(k).gen{indx},r];
                    end
                    unvisited = [];
                end
            end
        elseif object == 2
            while ~isempty(unvisited) 
                if ~isempty(unassigned) 
                    for l = 1:length(unassigned)
                        i = unassigned(l);
                        dists = vecnorm(X(unvisited,:)-X(centrals(i),:),2,2);
                        [~,indx] = min(dists);
                        r = unvisited(indx);
                        pre = pop(k).gen{i}(end);
                        pop(k).gen{i} = [pop(k).gen{i},r];
                        pop(k).ltours(i) = pop(k).ltours(i) + D(pre+1,r+1);
                        unvisited(indx) = [];
                        if isempty(unvisited) 
                            break;
                        end
                    end
                    toremove = [];
                    for l = 1:length(unassigned)
                        i = unassigned(l);
                        r = pop(k).gen{i}(end);
                        if pop(k).ltours(i)+D(r+1,1) + dmin > Dmax
                            toremove = [toremove, i];
                        end
                    end
                    unassigned = unassigned(~ismember(unassigned,toremove));
                else
                    % length(unassigned) == 0 & length(unvisited) > 0
                    len_unvisited = length(unvisited);
                    for i = 1:len_unvisited
                        r = unvisited(i);
                        indx = randi(m);
                        pop(k).gen{indx} = [pop(k).gen{indx},r];
                    end
                    unvisited = [];
                end
            end
        end
        pop(k).totalDist = 0;
        pop(k).feasibility = ones(m,1);
        pop(k).routes = cell(m,1);
        pop(k).numarcs = zeros(m,1);
        pop(k) = groupGAcost(pop(k),m,Dmax,object,regions,centralPs,InterRegionsPoints,NoofRegions,1,toursearch);
        dists = vertcat(pop(1:k-1).totalDist);
        if ~ismember(pop(k).totalDist, dists)
            k = k+1;   
        end
    end
end

end