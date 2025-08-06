function [Q1_K_loc, rhs_loc] = q1_local_matrices_rhs(P_q1_k, x, y, nx, ny, gaussian_nodes, gaussian_weigths, number_nodes)

% P_k is the local (at the k element) information of the connectivity matrix P
% x and y are the vector containing the x and the y coordinates of nodes of the whole grid
% nx and ny are the number of points on the x and y axis (in total) respectively
% hx and hy are the steplenght in the x and y direction respectively

    Q1_K_loc = zeros(4,4);
    
    rhs_loc = zeros(4,1);
    
    % the position f2(i,j) refers to the function evaluation
    % basis_finction_i at the j-th Gaussian point for the quadrature 
    f1 = zeros(4, number_nodes);
    g1 = zeros(4, number_nodes);
    
    f2 = zeros(4,number_nodes);
    g2 = zeros(4,number_nodes);
    
    force_func = zeros(number_nodes, number_nodes);
    
    position_first_node = P_q1_k(1,1);
    
    j = rem(position_first_node, ny);
    i = (position_first_node - j) / ny + 1;
    
    xi = x(i);
    yj = y(j);
    
    % meshsizes
    hy = y(j+1) - y(j);
	hx = x(i+1) - x(i);
    
    % evaluation of basis functions and their gradients
    for k = 1 : number_nodes
        % the functions are evaluated at each point on the x and y axis
        
        % Q1 basis functions
        [f1(1,k), g1(1,k)] = q1_basis_function_1(xi, yj, hx, hy, xi+hx/2+(hx/2)*gaussian_nodes(k), yj+hy/2+(hy/2)*gaussian_nodes(k));
        [f1(2,k), g1(2,k)] = q1_basis_function_2(xi, yj, hx, hy, xi+hx/2+(hx/2)*gaussian_nodes(k), yj+hy/2+(hy/2)*gaussian_nodes(k));
        [f1(3,k), g1(3,k)] = q1_basis_function_3(xi, yj, hx, hy, xi+hx/2+(hx/2)*gaussian_nodes(k), yj+hy/2+(hy/2)*gaussian_nodes(k));
        [f1(4,k), g1(4,k)] = q1_basis_function_4(xi, yj, hx, hy, xi+hx/2+(hx/2)*gaussian_nodes(k), yj+hy/2+(hy/2)*gaussian_nodes(k));
        
        
        % gradient of Q1 basis functions
        [f2(1,k), g2(1,k)] = q1_grad_basis_function_1(xi, yj, hx, hy, xi+hx/2+(hx/2)*gaussian_nodes(k), yj+hy/2+(hy/2)*gaussian_nodes(k));
        [f2(2,k), g2(2,k)] = q1_grad_basis_function_2(xi, yj, hx, hy, xi+hx/2+(hx/2)*gaussian_nodes(k), yj+hy/2+(hy/2)*gaussian_nodes(k));
        [f2(3,k), g2(3,k)] = q1_grad_basis_function_3(xi, yj, hx, hy, xi+hx/2+(hx/2)*gaussian_nodes(k), yj+hy/2+(hy/2)*gaussian_nodes(k));
        [f2(4,k), g2(4,k)] = q1_grad_basis_function_4(xi, yj, hx, hy, xi+hx/2+(hx/2)*gaussian_nodes(k), yj+hy/2+(hy/2)*gaussian_nodes(k));
        
        
        for k2 = 1 : number_nodes
            % force function evaluation
            force_func(k,k2) = f(xi+hx/2+(hx/2)*gaussian_nodes(k), yj+hy/2+(hy/2)*gaussian_nodes(k2));
        end
        
    end
    
    
    for i = 1 : 4
        for j = 1 : 4
            I2_x = 0;
            I2_y = 0;
            
            for k = 1 : number_nodes
                I2_y = I2_y + gaussian_weigths(k) * f2(i,k) * f2(j,k);
                I2_x = I2_x + gaussian_weigths(k) * g2(i,k) * g2(j,k);
            end

            I2_x = 0.5 * hx * I2_x;
            I2_y = 0.5 * hy * I2_y;

            Q1_K_loc(i,j) = hx * I2_y + hy * I2_x;
        end
        
        for k = 1 : number_nodes
            for k2 = 1 : number_nodes
            	rhs_loc (i,1) = rhs_loc(i,1) + gaussian_weigths(k) * gaussian_weigths(k2) * force_func(k,k2) * f1(i, k) * g1(i, k2);
            end
        end
        
        rhs_loc(i,1) = 0.25 * hx * hy * rhs_loc(i,1);
        
    end

end