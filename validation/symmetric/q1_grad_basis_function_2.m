function [x_component_grad_basis_function2, y_component_grad_basis_function2] = q1_grad_basis_function_2(xi, yj, hx, hy, x, y)

    % x_grad_component_basis_function2: component of the x variable to the
    % gradient of the basis function 2
    % y_grad_component_basis_function2: component of the y variable to the
    % gradient of the basis function 2

    
    x_component_grad_basis_function2 = hx^-1 * (1 + yj/hy - y/hy);
    
    y_component_grad_basis_function2 = -hy^-1 * (-xi/hx + x/hx);

end