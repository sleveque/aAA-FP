function [x_component_basis_function1, y_component_basis_function1] = q1_basis_function_1(xi, yj, hx, hy, x, y)

    % x_component_basis_function1: contribute of the x variable to the
    % basis function 1
    % y_component_basis_function1: contribute of the y variable to the
    % basis function 1

    
    x_component_basis_function1 = 1 + xi/hx - x/hx;
    
    y_component_basis_function1 = 1 + yj/hy - y/hy;

end