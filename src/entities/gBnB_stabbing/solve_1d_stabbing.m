function [overlap_maxnum, overlap_point] = solve_1d_stabbing(x_min, x_max)
% x_min: starting points of line segments
% x_max: endding points of line segments

    n = length(x_min);
    x_min = sort(x_min);
    x_max = sort(x_max);

    overlap_num = 0;
    overlap_maxnum = 0;
    
    i = 1;
    j = 1;
    
    while i <= n && j <= n
        if x_min(i) < x_max(j)
            overlap_num = overlap_num + 1;
            if overlap_num > overlap_maxnum
                overlap_maxnum = overlap_num;
                overlap_point = (x_min(i) + x_max(j)) / 2;
            end
            i = i + 1;
        else
            overlap_num = overlap_num - 1;
            j = j + 1;
        end
    end
    
end