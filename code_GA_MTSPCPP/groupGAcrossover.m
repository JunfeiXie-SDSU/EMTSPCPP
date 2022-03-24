function [child1,child2] = groupGAcrossover(indx_parents,pop,NoofRegions,m,object)
parent1 = pop(indx_parents(1));
parent2 = pop(indx_parents(2));
child1.gen = cell(m,1);
child2.gen = cell(m,1);

method = 1; % method for object 1

% pick the best trip from one parent
% if object == 1 % pick the feasible trip with the smallest ratio of length to the number of regions.
% otherwise, pick the shortest trip
%% generate child1
if object == 1
    if method == 1
        fi_indx = find(parent1.feasibility == 1);
        ratios = parent1.ltours./parent1.numarcs;
        if ~isempty(fi_indx)
            [~,inx] = min(ratios(fi_indx));
            child1.gen{1} = parent1.gen{fi_indx(inx)};
            visited = parent1.gen{fi_indx(inx)};
        else
            [~,inx] = min(ratios);
            child1.gen{1} = parent1.gen{inx};
            visited = parent1.gen{inx};
        end
    else
        fi_indx = find(parent1.feasibility == 1);
        if ~isempty(fi_indx)
            [~,inx] = max(parent1.ltours(fi_indx));
            child1.gen{1} = parent1.gen{fi_indx(inx)};
            visited = parent1.gen{fi_indx(inx)};
        else
            [~,inx] = min(parent1.ltours);
            child1.gen{1} = parent1.gen{inx};
            visited = parent1.gen{inx};
        end
    end
else
    [~,inx] = min(parent1.ltours);
    child1.gen{1} = parent1.gen{inx};
    visited = parent1.gen{inx};
end
unvisited = 1:NoofRegions;
unvisited = unvisited(~ismember(unvisited,visited));
if m == 2
    child1.gen{2} = unvisited;
elseif m > 2
    numarcs = zeros(m,1);
    sets = cell(m,1);
    for i = 1:m
        sets{i} = parent2.gen{i}(~ismember(parent2.gen{i},visited));
        numarcs(i) = length(sets{i});
    end
    [ix,~] = find(numarcs ~= 0);
    if length(ix) == m
        [~,li_inx] = min(numarcs);
        ix_temp = ix(~ismember(ix,li_inx));
        unassigned = sets{li_inx};
        while ~isempty(unassigned)
            [~,li_temp] = min(numarcs(ix_temp));
            r = randi(length(unassigned));
            sets{ix_temp(li_temp)} = [sets{ix_temp(li_temp)}, unassigned(r)];
            numarcs(ix_temp(li_temp)) = numarcs(ix_temp(li_temp)) + 1;
            unassigned(r) = [];
        end
        for j = 2:m
            child1.gen{j} = sets{ix_temp(j-1)};
        end
    else
        numzeros = m - length(ix);
        for j = 2:length(ix)+1
            child1.gen{j} = sets{ix(j-1)};
        end
        if numzeros > 1
            for j = length(ix)+2:m
                [lm, l_inx] = max(numarcs(ix));
                r = randi(lm);
                child1.gen{j} = sets{ix(l_inx)}(r);
                sets{ix(l_inx)}(r)=[];
                numarcs(ix(l_inx)) = numarcs(ix(l_inx))-1;
            end
            for j = 2:length(ix)+1
                child1.gen{j} = sets{ix(j-1)};
            end
        end
    end
end
    
%% generate child2
if object == 1
    if method == 1
        fi_indx = find(parent2.feasibility == 1);
        ratios = parent2.ltours./parent2.numarcs;
        if ~isempty(fi_indx)
            [~,inx] = min(ratios(fi_indx));
            child2.gen{1} = parent2.gen{fi_indx(inx)};
            visited = parent2.gen{fi_indx(inx)};
        else
            [~,inx] = min(ratios);
            child2.gen{1} = parent2.gen{inx};
            visited = parent2.gen{inx};
        end
    else
        fi_indx = find(parent2.feasibility == 1);
        if ~isempty(fi_indx)
            [~,inx] = max(parent2.ltours(fi_indx));
            child2.gen{1} = parent2.gen{fi_indx(inx)};
            visited = parent2.gen{fi_indx(inx)};
        else
            [~,inx] = min(parent2.ltours);
            child2.gen{1} = parent2.gen{inx};
            visited = parent2.gen{inx};
        end
    end
else
    [~,inx] = min(parent2.ltours);
    child2.gen{1} = parent2.gen{inx};
    visited = parent2.gen{inx};
end
unvisited = 1:NoofRegions;
unvisited = unvisited(~ismember(unvisited,visited));
if m == 2
    child2.gen{2} = unvisited;
elseif m > 2
    numarcs = zeros(m,1);
    sets = cell(m,1);
    for i = 1:m
        sets{i} = parent1.gen{i}(~ismember(parent1.gen{i},visited));
        numarcs(i) = length(sets{i});
    end
    [ix,~] = find(numarcs ~= 0);
    if length(ix) == m
        [~,li_inx] = min(numarcs);
        ix_temp = ix(~ismember(ix,li_inx));
        unassigned = sets{li_inx};
        while ~isempty(unassigned)
            [~,li_temp] = min(numarcs(ix_temp));
            r = randi(length(unassigned));
            sets{ix_temp(li_temp)} = [sets{ix_temp(li_temp)}, unassigned(r)];
            numarcs(ix_temp(li_temp)) = numarcs(ix_temp(li_temp)) + 1;
            unassigned(r) = [];
        end
        for j = 2:m
            child2.gen{j} = sets{ix_temp(j-1)};
        end
    else
        numzeros = m - length(ix);
        for j = 2:length(ix)+1
            child2.gen{j} = sets{ix(j-1)};
        end
        if numzeros > 1
            for j = length(ix)+2:m
                [lm, l_inx] = max(numarcs(ix));
                if lm == 0
                    lm = 0;
                end
                r = randi(lm);
                child2.gen{j} = sets{ix(l_inx)}(r);
                sets{ix(l_inx)}(r)=[];
                numarcs(ix(l_inx)) = numarcs(ix(l_inx))-1;
            end
            for j = 2:length(ix)+1
                child2.gen{j} = sets{ix(j-1)};
            end
        end
    end
end

end