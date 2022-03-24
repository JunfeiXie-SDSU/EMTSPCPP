function [region_order, subDistance, OptRoute] = findSubtour2(indx, regions, ...
    InterRegionsPoints,centralPs,NoofRegions,method)
global UAS_initP
% method = 'opt'; % use Fast-2opt to find the best tour
% method = 'nn'; % use nearest neighbor to determine region visiting order
% method = 'none';
N = length(indx);
    if N == 1
         region_order = [1, indx+1, 1];
        [OptRoute, subDistance] = FindTourWithRegionOrder2StepM(regions, region_order, InterRegionsPoints,centralPs);
    else
        subgroup = cell(N,1);
        subIRP = cell(N,1);
        for ii = 1:N
            subgroup{ii} = regions{indx(ii)};
            subIRP{ii} = InterRegionsPoints{indx(ii)};
        end
        subCentralPs = centralPs(indx,:);
        switch method
            case 'nn'  %nearest neighbor
            grids = [UAS_initP; subCentralPs];
            indices = [1,indx + ones(1,N)];
            [OptimalT,~] = tsp_nn(grids);
            [~, c] = find(OptimalT == 1);
            region_order1 = indices([OptimalT(c:end), OptimalT(1:c-1),1]);
            %% find the TSP-CPP tour
            [OptRoute1, subDistance1] = FindTourWithRegionOrder2StepM(regions, region_order1, InterRegionsPoints,centralPs);
            region_order2 = region_order1(end:-1:1);
            [OptRoute2, subDistance2] = FindTourWithRegionOrder2StepM(regions, region_order2, InterRegionsPoints,centralPs);
            if subDistance1 < subDistance2
                region_order = region_order1;
                subDistance = subDistance1;
                OptRoute = OptRoute1;
            else
                region_order = region_order2;
                subDistance = subDistance2;
                OptRoute = OptRoute2;
            end
            case 'opt'
            [OptRoute, subDistance,order_indx] = FastNNOpt(N, subgroup, subIRP,subCentralPs);
            region_order = [1, indx(order_indx)+ones(1,length(indx)),1];
            case 'none'
            region_order = [1, indx+ones(1,N), 1];
            [OptRoute, subDistance] = FindTourWithRegionOrder2StepM(regions, region_order, InterRegionsPoints,centralPs);
        end
    end
end
    