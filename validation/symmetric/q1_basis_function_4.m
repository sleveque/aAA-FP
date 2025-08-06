function [x_component_basis_function4, y_component_basis_function4] = q1_basis_function_4(xi, yj, hx, hy, x, y)

    % x_component_basis_function4: contribute of the x variable to the
    % basis function 4
    % y_component_basis_function4: contribute of the y variable to the
    % basis function 4

    
    x_component_basis_function4 = 1 + xi/hx - x/hx;
    
    y_component_basis_function4 = -yj/hy + y/hy;

end