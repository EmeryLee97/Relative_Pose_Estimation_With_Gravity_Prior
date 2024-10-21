function [t_opt_3d, theta_opt] = solve_BnB_stabbing(pts_1, pts_2, R_v, epsilon)
% pts_1: matching points in the 1st image
% pts_2: matching points in the 2nd image
% R_v: rotation axis
% epsilon: tolerance
    
%     ss=get(0,'screensize');
%     
%     len_x=0.45*ss(3);
%     len_y=0.45*ss(4);
%     
%     figure;
%     set(gcf,'position',[1 1 len_x len_y]);
%     title('Evolution of the bounds of exp-map')
%     hold on
%     grid
%     h_l=animatedline('color','r');
%     h_u=animatedline('color','b');
%     legend('upper bound','lower bound')
% 
%     iter_times = 0;
    
    % Centered at (0, 0), radius = pi/2
    solution_domain = [0; 0; pi/2];
    branch_assembled = [];
    
    best_branch = solution_domain;

    global_lower_bound = 0;
    global_upper_bound = size(pts_1, 2);

%     % The rotation that maps g1 to g2 with min geodesic motion
%     angle_g1_g2 = acos(g_1' * g_2);
%     axis_g1_g2 = cross(g_1, g_2);
%     R_g1_g2 = axang2rotm([axis_g1_g2', angle_g1_g2]);
% 
%     g2_cross = [0 -g_2(3) g_2(2); g_2(3) 0 -g_2(1); -g_2(2) g_2(1) 0];
% 
%     % a, b, c are (3, N) constant matricies
%     a = cross(pts_2, ((eye(3) + g2_cross^2) * R_g1_g2 * pts_1));
%     b = cross(pts_2, g2_cross * R_g1_g2 * pts_1);
%     c = -cross(pts_2, g2_cross^2 * R_g1_g2 * pts_1);

    K=[0 -R_v(3) R_v(2);
        R_v(3) 0 -R_v(1);
        -R_v(2) R_v(1) 0]';
    
    b = cross(pts_2, K * pts_1);
    c = -cross(pts_2, K^2 * pts_1);
    a = cross(pts_2, pts_1) + cross(pts_2, K^2 * pts_1);


    while global_upper_bound > global_lower_bound

        % Updated properties during each iteration
        new_radius = best_branch(end) / 2;
        new_upper_bounds = zeros(1, 4);
        new_lower_bounds = zeros(1, 4);
        theta_values = zeros(1, 4);

        % Center points of the new branches
        cp_1 = [best_branch(1) + new_radius; best_branch(2) + new_radius];
        cp_2 = [best_branch(1) - new_radius; best_branch(2) - new_radius];
        cp_3 = [best_branch(1) - new_radius; best_branch(2) + new_radius];
        cp_4 = [best_branch(1) + new_radius; best_branch(2) - new_radius];
        current_branch = [cp_1, cp_2, cp_3, cp_4; new_radius * ones(1, 4)]; % 3*4 matrix

        for i = 1:4
            [new_upper_bounds(i), new_lower_bounds(i), theta_values(i)] = ...
                get_bound_theta917(a, b, c, current_branch(:, i), epsilon);
        end
        

        % Assemble the branches with calculated bounds, 5*k matrix
        % [center_1, center_2, half_side_length, lower_bound, upperbound]
        branch_assembled = [branch_assembled, [current_branch; new_lower_bounds; new_upper_bounds]];

        % Update the optimal values when a bigger lower bound appears
        [new_lower_bound, ind_lower] = max(new_lower_bounds);
        if global_lower_bound < new_lower_bound
            global_lower_bound = new_lower_bound;
            t_opt_2d = current_branch(1:2, ind_lower);
            theta_opt = theta_values(ind_lower);
        end

        % Let the maximum upper bound be the global upper bound
        [global_upper_bound, ind_global_upper] = max(branch_assembled(end, :));
        %dbstop in solve_BnB_stabbing at 97 if global_upper_bound<50
        %dbstop in get_bound_theta917 at 148 if global_upper_bound<50
        
        % Go deep into the best branch which contains the maximum upper bound
        best_branch = branch_assembled(1:3, ind_global_upper);
        branch_assembled(:, ind_global_upper) = [];

        % Prune branches whose upper bound is less than the global lower bound
        branch_assembled(:, branch_assembled(end, :) < global_lower_bound) = [];

%         iter_times = iter_times + 1;
% 
%         addpoints(h_l,iter_times,global_upper_bound);
%         addpoints(h_u,iter_times,global_lower_bound);
%         drawnow    

    end
    
    % convert t_2d back to t_3d
    t_opt_3d = exponential_mapping_2dto3d(t_opt_2d);

end