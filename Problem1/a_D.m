% Guard set
function [inside] = a_D(x) 
    if x(2) >= 1
        inside = 1;
    else
        inside = 0;
    end
end
