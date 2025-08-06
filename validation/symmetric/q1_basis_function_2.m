function [x_component_basis_function2, y_component_basis_function2] = q1_basis_function_2(xi, yj, hx, hy, x, y)

    % x_component_basis_function2: contribute of the x variable to the
    % basis function 2
    % y_component_basis_function2: contribute of the y variable to the
    % basis function 2
    

    x_component_basis_function2 = -xi/hx + x/hx;
    
    y_component_basis_function2 = 1 + yj/hy - y/hy;

end