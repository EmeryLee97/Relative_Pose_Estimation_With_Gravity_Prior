function [upper_bound, lower_bound, theta] = get_bound_theta905(a, b, c, branch_2d, epsilon)

    % This function finds the most overlapping theta
    % a, b, c: (3, n) matrices derived from paper
    % branch: [branch_center; branch_half_length]
    % epsilon: tolerence

    % Branch center & branch half length
    t_c = exponential_mapping_2dto3d(branch_2d(1:2));
    delta = branch_2d(3);

    % |t'(a+b*sin(theta)+c*cos(theta))| <= epsilon
    % |A+B*sin(theta)+C*cos(theta)| <= epsilon
    A = t_c' * a; % (1, n)
    B = t_c' * b; % (1, n)
    C = t_c' * c; % (1, n)

    % |cos(theta-phi)+k1| <= epsilon_tilde
    % with cos(phi) = B / sqrt(A^2+B^2)
    denominator = sqrt(B.^2 + C.^2); % (1, n)
    k1 = A ./ denominator; % (1, n)
    epsilon_tilde = epsilon ./ denominator; % (1, n)
    phi = atan(B ./ C); % (1, n)

    % case 1: |theta-phi| <= pi
    % k2 <= cos(theta-phi) <= k3
    k2 = max(-1, -epsilon_tilde-k1); % (1, n)
    k3 = min(1, epsilon_tilde-k1); % (1, n)

    % theta ~ [phi-arccos(k2), phi-arccos(k3)] || 
    %         [phi+arccos(k3), phi+arccos(k2)] && 
    %         [-pi, pi]
    k4 = phi - acos(k2); k5 = phi - acos(k3); % (1, n)
    k6 = phi + acos(k3); k7 = phi + acos(k2); % (1, n)

    % case 2: pi <= |theta-phi| <= 2*pi
    % k2 <= cos(2*pi-(theta-phi)) <= k3

    % theta ~ [-2*pi+phi+arccos(k3), -2*pi+phi+arccos(k2)] ||
    %         [2*pi+phi-arccos(k2), 2*pi+phi-arccos(k3)] &&
    %         [-pi, pi]

    theta_set = [k4, k6, -2*pi+k6, 2*pi+k4; 
                 k5, k7, -2*pi+k7, 2*pi+k5];

    theta_set(:, theta_set(1, :) > pi) = [];
    theta_set(:, theta_set(2, :) < -pi) = [];
    theta_set(1, theta_set(1, :) < -pi) = -pi;
    theta_set(2, theta_set(2, :) > pi) = pi; % (2, m)

    % interval stabbing
    [lower_bound, theta] = solve_1d_stabbing(theta_set(1, :), theta_set(2, :), "middle");
    %[lower_bound, theta] = Max_Stabbing(theta_set(1, :), theta_set(2, :));

    % k7 = a+b*sin(theta)+c*cos(theta)
    k7 = a + b .* sin(theta) + c .* cos(theta); % (3, n)

%     % k8 = angle(t_c, k7)
%     k8 = acos((t_c' * k7) ./ (vecnorm(t_c) .* vecnorm(k7))); % (1, n), [0, pi]
%     
%     % k9 = angle(t, k7), with cos(k9) = max
%     k9 = k8;
%     k9(k8 - sqrt(2) * delta < 0) = 0; % (1, n)
%     k9(k8 - sqrt(2) * delta >= 0) = k8(k8 - sqrt(2) * delta >= 0) - sqrt(2) * delta;
% 
%     % k10 = angle(t, k7), with cos(k10) = min
%     k10 = k8;
%     k10(k8 + sqrt(2) * delta > pi) = pi; % (1, n)
%     k10(k8 + sqrt(2) * delta <= pi) = k8(k8 + sqrt(2) * delta <= pi) + sqrt(2) * delta;
% 
%     % k11 = d_max
%     k11 = vecnorm(k7) .* vecnorm(t_c) .* cos(k9); % (1, n)
%     
%     % k12 = d_min
%     k12 = vecnorm(k7) .* vecnorm(t_c) .* cos(k10); % (1, n)
% 
%     % upper_bound
%     upper_bound = sum((k11 >= -epsilon) .* (k12 <= epsilon));
end