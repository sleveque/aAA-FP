function v_boundary_K = construct_boundary_conditions_for_v(a, b, c, d, nx_exp, ny_exp, sigma_x, sigma_y, K_boundary, FE_order)

    if FE_order == 2
        % numbers of intervals on the x and y axis
        n = 2 ^ nx_exp + 1;
        m = 2 ^ ny_exp + 1;
    elseif FE_order == 3
        % numbers of intervals on the x and y axis
        n = 2 ^ (nx_exp + 1) + 1;
        m = 2 ^ (ny_exp + 1) + 1;
    end
    
    % number of (interior) points on the x and y axis
    nx = n - 2;
	ny = m - 2;
    
    
    % construction of the vectors containing the coordinates of the grids
    if sigma_x == 0 && sigma_y == 0
        x_q1 = linspace(a, b, n);
        y_q1 = linspace(c, d, m);
    else
        x_q1 = zeros(n,1);
        y_q1 = zeros(m,1);
        
        x_q1(1:(n-1)/2+1) = linspace(a,sigma_x,(n+1)/2);
        x_q1((n-1)/2+1:end) = linspace(sigma_x,b,(n+1)/2); % overwrite middle point
        y_q1(1:(m-1)/2+1) = linspace(c,sigma_y,(m+1)/2);
        y_q1((m-1)/2+1:end) = linspace(sigma_y,d,(m+1)/2); % overwrite middle point
    end

    
    % construction of the vectors containing the boundary conditions
    left_boundary = zeros(m,1);
    low_boundary = zeros(nx,1);
    right_boundary = zeros(m,1);
    up_boundary = zeros(nx,1);
    
	for i = 1 : m
        left_boundary(i,1) = bound_value_left(x_q1(1),y_q1(i));
        right_boundary(i,1) = bound_value_right(x_q1(m),y_q1(i));
	end
        
	for j = 2 : n-1
    	low_boundary(j-1,1) = bound_value_down(x_q1(j),y_q1(1));
        up_boundary(j-1,1) = bound_value_up(x_q1(j),y_q1(m));
	end

    boundary = [left_boundary; low_boundary; up_boundary; right_boundary];
    
    v_boundary_K = K_boundary * boundary;

end