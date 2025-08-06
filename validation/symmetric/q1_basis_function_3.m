function [x_component_basis_function3, y_component_basis_function3] = q1_basis_function_3(xi, yj, hx, hy, x, y)

    % x_component_basis_function3: contribute of the x variable to the
    % basis function 3
    % y_component_basis_function3: contribute of the y variable to the
    % basis function 3

    
    x_component_basis_function3 = -xi/hx + x/hx;
    
    y_component_basis_function3 = -yj/hy + y/hy;

end