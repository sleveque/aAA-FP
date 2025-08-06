function [x_component_grad_basis_function3, y_component_grad_basis_function3] = q1_grad_basis_function_3(xi, yj, hx, hy, x, y)

    % x_grad_component_basis_function3: component of the x variable to the
    % gradient of the basis function 3
    % y_grad_component_basis_function3: component of the y variable to the
    % gradient of the basis function 3

    
    x_component_grad_basis_function3 = hx^-1 * (-yj/hy + y/hy);
    
    y_component_grad_basis_function3 = hy^-1 * (-xi/hx + x/hx);

end