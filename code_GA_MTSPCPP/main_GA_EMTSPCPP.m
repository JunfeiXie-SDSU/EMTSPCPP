%=========================================================================
% Genetic Algorithm for EMTSP-CPP (Energy Constrainted Multiple TSP-CPP)
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
load 'FiftyRegionsCase.mat'  % load the map of the regions
global UAS_range 
global UAS_initP
UAS_range = [1.5,3];  % sensing range of the UAVs;
UAS_initP = [100,100]; %initial position of the UAVs. All UAVs depart from the same location;
m = 6; % number of UAVs
% if object=1, minimize the total length; if object = 2, minimize the
% longest tour; 
object = 1; 
Dmax = 1000; % maximum distance a UAV can travel
tic;
%% find candidate coverage paths for each region
InterRegionsPoints = FindAllBFPPaths3(NoofRegions, regions);

%% find the central point for each region
centralPs = zeros(NoofRegions, 2);
for i = 1:NoofRegions
    rect = regions{i};
    centralPs(i,:) = mean(rect,1);
end
time_load = toc;

%% find the TSP-CPP tour
popSize = 100;  % population size
numIter = 1e8; % maximum number of iterations to run
pm = 0.05; % mutation probability
TermCrit = 10;  % algorithm is considered to converge if solution found does not change for TermCrit consecutive iterations 
Er = 0.1;  % elitism rate

% term_method = 'count';
term_method = 'time';
time_max = 60 - time_load;

%% initialization
toursearchi = 'none'; %% no algorithm applied
%% group GA
toursearch = 'opt'; %% Fast2OPT

tic;
%% initialization
pop_init = InitGroup(m, NoofRegions, popSize, Dmax,regions,InterRegionsPoints,centralPs,object,toursearchi);
time_init = toc;
dists = vertcat(pop_init.totalDist);

time_max = time_max - time_init;

tic;
%% Apply group GA
[globalMin_group,OptRoute_group,iter_group] = groupGA(pop_init,m,NoofRegions,object,Dmax,regions,...
    centralPs,InterRegionsPoints,numIter,popSize,Er,pm,TermCrit,term_method,time_max,toursearch);
time_group = toc;


%% plot figure
colors =  ['b','k','y','c','g','r','m'];
% %%%%% plot the trajectory for group GA
figure;
hold on;
for i = 1:NoofRegions
    rectangle = regions{i};
    fill(rectangle(:,1),rectangle(:,2),[0.9,0.9,0.9])
end

ltours_group = zeros(m,1);
for j = 1:m
    for i = 1:length(OptRoute_group{j}(:,1)) -1
        line(OptRoute_group{j}(i:i+1,1), OptRoute_group{j}(i:i+1,2), 'Color', colors(j), 'Marker','.', 'LineWidth', 1.5);
        ltours_group(j) = ltours_group(j) + norm(OptRoute_group{j}(i,:)-OptRoute_group{j}(i+1,:));
    end
end
ltours_group
