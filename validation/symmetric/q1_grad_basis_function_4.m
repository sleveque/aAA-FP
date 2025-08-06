function [x_component_grad_basis_function4, y_component_grad_basis_function4] = q1_grad_basis_function_4(xi, yj, hx, hy, x, y)

    % x_grad_component_basis_function4: component of the x variable to the
    % gradient of the basis function 4
    % y_grad_component_basis_function4: component of the y variable to the
    % gradient of the basis function 4

    
    x_component_grad_basis_function4 = -hx^-1 * (-yj/hy + y/hy);
    
    y_component_grad_basis_function4 = hy^-1 * (1 + xi/hx - x/hx);

end