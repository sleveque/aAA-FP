function [x_component_grad_basis_function1, y_component_grad_basis_function1] = q1_grad_basis_function_1(xi, yj, hx, hy, x, y)

    % x_grad_component_basis_function1: component of the x variable to the
    % gradient of the basis function 1
    % y_grad_component_basis_function1: component of the y variable to the
    % gradient of the basis function 1

    
    x_component_grad_basis_function1 = -hx^-1 * (1 + yj/hy - y/hy);
    
    y_component_grad_basis_function1 = -hy^-1 * (1 + xi/hx - x/hx);

end