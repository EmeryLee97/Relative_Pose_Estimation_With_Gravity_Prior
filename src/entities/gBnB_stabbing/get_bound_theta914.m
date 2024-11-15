function [upper_bound, lower_bound, theta] = get_bound_theta914( ...
    a, b, c, branch_2d, epsilon)

    % Branch center & branch half length
    t_c = exponential_mapping_2dto3d(branch_2d(1:2));
    delta = branch_2d(3);

    % sample theta between [-sqrt(2)*delta, sqrt(2)*delta]
    % sample step = 0.02 (0.02rad = 1.14deg)
    gamma = -sqrt(2)*delta : 0.001 : sqrt(2)*delta; % (1, m)

    % norm of vectors a, b, c
    a_norm = vecnorm(a); % (1, n)
    b_norm = vecnorm(b); % (1, n)
    c_norm = vecnorm(c); % (1, n)

    % angles between tc and a, b, c
    angle_a_tc = acos(t_c' * a ./ a_norm); % (1, n)
    angle_b_tc = acos(t_c' * b ./ b_norm); % (1, n)
    angle_c_tc = acos(t_c' * c ./ c_norm); % (1, n)

    % angles between t and a, b, c               % a_gamma:  a(1) a(2) a(3) ... a(n)
    % angle_t_tc <= sqrt(2) * delta              % gamma(1): 
    % m = length(gamma), n = size(a, 2);         % gamma(2):
    angle_a_t = angle_a_tc + gamma'; % (m, n)    % gamma(3):     a+gamma
    angle_b_t = angle_b_tc + gamma'; % (m, n)    %  ...... 
    angle_c_t = angle_c_tc + gamma'; % (m, n)    % gamma(m):

    % A, B, C are cross products of a, b, c and t
    A = a_norm .* cos(angle_a_t); % (m, n)
    B = b_norm .* cos(angle_b_t); % (m, n)
    C = c_norm .* cos(angle_c_t); % (m, n)

    denominator = sqrt(B.^2 + C.^2); % (m, n)

    % sin(phi) > 0, cos(phi) > 0: (0, pi/2) 第一象限
    % sin(phi) > 0, cos(phi) < 0: (pi/2, pi) 第二象限
    % sin(phi) < 0, cos(phi) < 0: (-pi, -pi/2) 第三象限
    % sin(phi) < 0, cos(phi) > 0: (-pi/2, 0) 第四象限
    sin_phi = B ./ denominator;
    cos_phi = C ./ denominator;
    phi = atan(B ./ C); % [-1/2*pi, 1/2*pi]
    for i = 1:size(phi, 1)
        for j = 1:size(phi, 2)
            if phi(i, j) >= 0
                if sin_phi(i, j) < 0 && cos_phi(i, j) < 0
                    phi(i, j) = phi(i, j) - pi;
                end
            else
                if sin_phi(i, j) > 0 && cos_phi(i, j) < 0
                    phi(i, j) = phi(i, j) + pi;
                end
            end
        end
    end

    A_tilde = A ./ denominator; % (m, n)
    epsilon_tilde = epsilon ./ denominator; % (m, n)

    % k1 <= cos(theta - phi) <= k2, note that k2 might < -pi
    k1 = max(-1, -epsilon_tilde - A_tilde); % (m, n)
    k2 = min(1, epsilon_tilde - A_tilde); % (m, n)

    k2(k2 < -1) = -1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k1(k1 > 1) = 1;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    acos_k1 = acos(k1); % (m, n)
    acos_k2 = acos(k2); % (m, n)

    % case 1: |theta - phi| <= pi
    % theta ~ [phi - arccos(k1), phi - arccos(k2)] ||
    %         [phi + arccos(k2), phi + arccos(k1)] && [-pi, pi]

    theta_min_1 = min(phi - acos_k1, [], 1); % (1, n)
    theta_max_1 = max(phi - acos_k2, [], 1); % (1, n)

    theta_min_2 = min(phi + acos_k2, [], 1); % (1, n)
    theta_max_2 = max(phi + acos_k1, [], 1); % (1, n)

    % theta_max_1不会大于pi，但会小于-pi，
    % theta_min_2不会小于-pi，但会大于pi
%     if sum(theta_max_1 > pi) ~= 0 || sum(theta_min_2 < -pi) ~= 0
%         disp("error");
%     end

%     if sum(theta_max_1 < -pi) ~= 0  
%         disp("error1");
%     end
% 
%     if sum(theta_min_2 > pi) ~= 0
%         disp("error2");
%     end

%       if sum(theta_min_1 <= -2*pi) ~= 0
%           disp("error1");
%       end
% 
%       if sum(theta_max_2 >= 2*pi) ~= 0
%           disp("error2");
%       end

    if sum(theta_min_1+2*pi < theta_max_1) > 0 
        disp("error")
    end

    % case 2: pi <= |theta - phi| <= 2*pi
    % theta ~ [phi + arccos(k2) - 2*pi, phi + arccos(k1) - 2*pi] ||
    %         [phi - arccos(k1) + 2*pi, phi - arccos(k2) + 2*pi] && [-pi, pi]
    
%     theta_min_3 = min(phi + acos_k2, [], 1) - 2*pi;
%     theta_max_3 = max(phi + acos_k1, [], 1) - 2*pi;
% 
%     theta_min_4 = min(phi - acos_k1, [], 1) + 2*pi;
%     theta_max_4 = max(phi - acos_k2, [], 1) + 2*pi;
% 
%     idx_2 = theta_min_4 < theta_max_3;
%     theta_min_4(idx_2) = theta_max_3(idx_2) + 1e-8;

    % merge intervals
    % 第一层if讨论case1, 第二层if讨论case2
    theta_set = [];
    for i = 1 : length(theta_min_1)
        if theta_max_1(i) >= theta_min_2(i) % 两个区间相交
            if theta_min_1(i) >= -pi && theta_max_2(i) <= pi % 左右都不超
                theta_set = [theta_set, [theta_min_1(i); theta_max_2(i)]];
            elseif theta_min_1(i) < -pi && theta_max_2(i) <= pi % 左超右不超
                left = max(theta_min_1(i)+2*pi, theta_max_2(i));
                theta_set = [theta_set, [-pi, left; theta_max_2(i), pi]];
            elseif theta_min_1(i) >= -pi && theta_max_2(i) > pi % 右超左不超
                right = min(theta_min_1(i), theta_max_2(i)-2*pi);
                theta_set = [theta_set, [-pi, theta_min_1(i); right, pi]];
            else % 左右都超
                theta_set = [theta_set, [-pi; pi]];
            end
        else % 两个区间不相交
            if theta_min_1(i) >= -pi && theta_max_2(i) <= pi % 左右都不超
                theta_set = [theta_set, [theta_min_1(i), theta_min_2(i); theta_max_1(i), theta_max_2(i)]];
            elseif theta_min_1(i) < -pi && theta_max_2(i) <= pi % 左超右不超
                if theta_max_1(i)+2*pi < theta_min_2(i) % 左区间+2pi与右区间不重叠
                    left = theta_min_1(i)+2*pi;
                    right = theta_max_1(i)+2*pi;
                    theta_set = [theta_set, [left, theta_min_2(i); right, theta_max_2(i)]];
                else % 左区间+2pi与右区间重叠
                    if theta_max_1(i) <= -pi
                        left = min(theta_min_1(i)+2*pi, theta_min_2(i));
                        right = max(theta_max_1(i)+2*pi, theta_max_2(i));
                        theta_set = [theta_set, [left; right]];
                    else
                        left = min(theta_min_1(i)+2*pi, theta_min_2(i)); %%%%%%%%这里不对
                        theta_set = [theta_set, [-pi, left; theta_max_1(i), pi]];
                    end
                end

                %theta_set = [theta_set, [-pi, left; theta_max_1(i), pi]];
            elseif theta_min_1(i) >= -pi && theta_max_2(i) > pi % 右超左不超
                if theta_min_2(i)-2*pi > theta_max_1(i) % 右区间-2pi与左区间不重叠
                    left = theta_min_2(i)-2*pi;
                    right = theta_max_2(i)-2*pi;
                    theta_set = [theta_set, [theta_min_1(i), left; theta_max_1(i), right]];
                else % 右区间-2pi与左区间重叠
                    if theta_min_2(i) >= pi
                        left = min(theta_min_1(i), theta_min_2(i)-2*pi);
                        right = max(theta_max_1(i), theta_max_2(i)-2*pi);
                        theta_set = [theta_set, [left; right]];
                    else
                        right = max(theta_max_2(i)-2*pi, theta_max_1(i));
                        theta_set = [theta_set, [-pi, theta_min_2(i); right, pi]];
                    end
                end
            else % 左右都超
                if theta_max_1(i) <= -pi && theta_min_2(i) >= pi
                    theta_set = [theta_set, [theta_min_1(i)+2*pi, theta_min_2(i)-2*pi; theta_max_1(i)+2*pi, theta_max_2(i)-2*pi]];
                elseif theta_max_1(i) <= -pi && theta_min_2(i) < pi 
                    if theta_max_1(i)+2*pi >= theta_min_2(i)
                        left = min(theta_min_1(i)+2*pi, theta_min_2(i));
                        theta_set = [theta_set, [-pi, left; theta_max_2(i)-2*pi, pi]];
                    else
                        theta_set = [theta_set, [-pi, theta_min_1(i)+2*pi, theta_min_2(i); theta_max_2(i)-2*pi, theta_max_1(i)+2*pi, pi]];
                    end
                elseif theta_max_1(i) > -pi && theta_min_2(i) >= pi
                    if theta_min_2(i)-2*pi <= theta_max_1(i)
                        right = max(theta_max_1(i), theta_max_2(i)-2*pi);
                        theta_set = [theta_set, [-pi, theta_min_1(i)+2*pi; right, pi]];
                    else
                        theta_set = [theta_set, [-pi, theta_min_2(i)-2*pi, theta_min_1(i)+2*pi; theta_max_1(i), theta_max_2(i)-2*pi, pi]];
                    end
                else
                    left = min(theta_min_1(i)+2*pi, theta_min_2(i));
                    right = max(theta_max_2(i)-2*pi, theta_max_1(i));
                    theta_set = [theta_set, [-pi, left; right, pi]];
                end
            end
        end
    end

    % get upper bound
    [upper_bound, theta] = solve_1d_stabbing(theta_set(1, :), theta_set(2, :), "middle");

    % get lower bound
    d = t_c' * (a + b .* sin(theta) + c .* cos(theta));
    lower_bound = sum(abs(d) <= epsilon);
end

