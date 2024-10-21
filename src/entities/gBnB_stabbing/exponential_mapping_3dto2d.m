function point_2d = exponential_mapping_3dto2d(point_3d)

    idx = point_3d(1, :) < 0;
    point_3d(:, idx) = -point_3d(:, idx);

    rho = abs(acos(point_3d(1, :)));
    
    dirc = point_3d(2:3, :) ./ sin(rho);
    point_2d = rho .* dirc;
    
end