function [child] = groupGAmutation(child,m,Pm)

if m > 1
if rand < Pm
    % randomly pick two trips
    t1 = randi(m);
    t2 = randi(m);
    while t1 == t2
        t2 = randi(m);
    end
    
    ltrip1 = length(child.gen{t1});
    ltrip2 = length(child.gen{t2});

    ix1 = randi(ltrip1);
    ix2 = randi(ltrip2);
    
    temp = child.gen{t2}(ix2);
    child.gen{t2}(ix2) = child.gen{t1}(ix1);
    child.gen{t1}(ix1) = temp;
    
end
end

end