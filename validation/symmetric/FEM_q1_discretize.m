function [K, K_boundary, rhs] = FEM_q1_discretize(a, b, c, d, nx_exp, ny_exp, number_nodes, sigma_x, sigma_y)

% Discretization for the Finite Element Metod using Q1-elements
%
% input
%
% [a,b]x[c,d]      domain of integration
% nx_exp           2^nx_exp is the number of intervals on the x-axis
% ny_exp           2^ny_exp is the number of intervals on the y-axis
% number_nodes     number of nodes for the Gaussian quadrature
%
% output
%
% K                stiffness matrix
% K_boundary       boundary conditions of discretization of differential operator
% rhs              right hand side


    % numbers of intervals on the x and y axis
    n = 2 ^ nx_exp + 1;
    m = 2 ^ ny_exp + 1;
    
    % number of (interior) points on the x and y axis
    nx = n - 2;
    ny = m - 2;
    
    % dimension of the (interior) grid
    dim = nx * ny;
    
    
    % construction of the vectors containing the coordinates of the grids
    if sigma_x == 0 && sigma_y == 0
        x_q1 = linspace(a, b, n);
        y_q1 = linspace(c, d, m);
    else
        x_q1 = zeros(n,1);
        y_q1 = zeros(m,1);
        
        x_q1(1:2^(nx_exp-1)+1) = linspace(a,sigma_x,(n+1)/2);
        x_q1(2^(nx_exp-1)+1:end) = linspace(sigma_x,b,(n+1)/2); % overwrite middle point
        y_q1(1:2^(ny_exp-1)+1) = linspace(c,sigma_y,(m+1)/2);
        y_q1(2^(ny_exp-1)+1:end) = linspace(sigma_y,d,(m+1)/2); % overwrite middle point
    end
    
    
    % construction of the connectivity matrices for Q1 elements
    P_q1 = q1_connectivity_matrix(n-1, m-1);

    
    % construction of the Gaussian nodes and weigths
    [gaussian_nodes, gaussian_weights] = gaussian_nodes_and_weights(number_nodes);
    
    
    % construction of the matrix and rhs
    K = sparse(dim, dim);
    
    rhs = zeros(dim, 1);
    
    
	% construction of the matrices for boundary conditions
    num_boundary = 2 * (nx + ny + 2);
    
    K_boundary = sparse(dim, num_boundary);
        
    
    % iteration on the left part of the tessellation
    P_q1_k= P_q1(1,:);
    
    [Q1_K_loc, rhs_1_loc] = q1_local_matrices_rhs(P_q1_k, x_q1, y_q1, n, m, gaussian_nodes, gaussian_weights, number_nodes);

    K = K + sparse(1,1,Q1_K_loc(3,3),dim,dim);
    
    rhs(1,1) = rhs(1,1) + rhs_1_loc(3,1);
    
    index_i_boundary = ones(3,1);
    index_j_boundary = [1;2;m+1];
    
    element_K_v_boundary = [Q1_K_loc(3,1),Q1_K_loc(3,4),Q1_K_loc(3,2)];
    
    K_boundary = K_boundary + sparse(index_i_boundary,index_j_boundary,element_K_v_boundary,dim,num_boundary);
    
    for k = 2 : m - 2
        P_q1_k= P_q1(k,:);
        
        index = zeros(1,2);
        index(1,1) = P_q1_k(1,2) - ny - 3;
        index(1,2) = P_q1_k(1,3) - ny - 3;
        
        [Q1_K_loc, rhs_1_loc] = q1_local_matrices_rhs(P_q1_k, x_q1, y_q1, n, m, gaussian_nodes, gaussian_weights, number_nodes);
        
        index_i = kron(index',ones(2,1));
        index_j = kron(ones(2,1),index');
        
        element_K_v1 = reshape(Q1_K_loc([2,3],[2,3])',4,1);
    
        K = K + sparse(index_i,index_j,element_K_v1,dim,dim);
        
        for i = 1 : 2
            rhs(index(1,i),1) = rhs(index(1,i),1) + rhs_1_loc(i+1,1);
        end

        index_i_boundary = kron(index',ones(2,1));
        index_j_boundary = kron(ones(2,1),[P_q1_k(1);P_q1_k(4)]);
        
        element_K_v_boundary = reshape(Q1_K_loc([2,3],[1,4])',4,1);

        K_boundary = K_boundary + sparse(index_i_boundary,index_j_boundary,element_K_v_boundary,dim,num_boundary);

    end
    
    
    P_q1_k= P_q1(m-1,:);
    
    [Q1_K_loc, rhs_1_loc] = q1_local_matrices_rhs(P_q1_k, x_q1, y_q1, n, m, gaussian_nodes, gaussian_weights, number_nodes);
    
    K = K + sparse(ny,ny,Q1_K_loc(2,2),dim,dim);
    
    rhs(ny,1) = rhs(ny,1) + rhs_1_loc(2,1);
    
    index_i_boundary = [ny;ny;ny];
    index_j_boundary = [m-1;m;m+nx+1];
        
    element_K_v_boundary = [Q1_K_loc(2,1),Q1_K_loc(2,4),Q1_K_loc(2,3)];

    K_boundary = K_boundary + sparse(index_i_boundary,index_j_boundary,element_K_v_boundary,dim,num_boundary);
    
	% iteration on the interior of the tessellation
    for k = 1 : n - 3
        P_q1_k= P_q1(k * (m - 1) + 1,:);
    
        [Q1_K_loc, rhs_1_loc] = q1_local_matrices_rhs(P_q1_k, x_q1, y_q1, n, m, gaussian_nodes, gaussian_weights, number_nodes);
        
        index = zeros(1,2);
        index(1) = P_q1_k(1,3) - ny - 3 - 2 * k;
        index(2) = P_q1_k(1,4) - ny - 1 - 2 * k;
        
        index_boundary = zeros(1,2);
        index_boundary(1) = ny+2+k;
        index_boundary(2) = ny+3+k;
        
        index_i = kron(index',ones(2,1));
        index_j = kron(ones(2,1),index');
        
        element_K_v1 = reshape(Q1_K_loc([3,4],[3,4])',4,1);
    
        K = K + sparse(index_i,index_j,element_K_v1,dim,dim);
        
        rhs(index(1),1) = rhs(index(1),1) + rhs_1_loc(3,1);
        rhs(index(2),1) = rhs(index(2),1) + rhs_1_loc(4,1);
        
        index_i_boundary = kron(index',ones(2,1));
        index_j_boundary = kron(ones(2,1),[index_boundary(1);index_boundary(2)]);
        
        element_K_v_boundary = reshape(Q1_K_loc([3,4],[1,2])',4,1);

        K_boundary = K_boundary + sparse(index_i_boundary,index_j_boundary,element_K_v_boundary,dim,num_boundary);
        
        for k2 = 2 : m -2
            P_q1_k= P_q1(k * (m - 1) + k2,:);
    
            [Q1_K_loc, rhs_1_loc] = q1_local_matrices_rhs(P_q1_k, x_q1, y_q1, n, m, gaussian_nodes, gaussian_weights, number_nodes);
        
            index = zeros(1,4);
            index(1) = P_q1_k(1,1) - ny - 1 - 2 * k;
            index(2) = P_q1_k(1,2) - ny - 3 - 2 * k;
            index(3) = P_q1_k(1,3) - ny - 3 - 2 * k;
            index(4) = P_q1_k(1,4) - ny - 1 - 2 * k;
            
            index_i = kron(index',ones(4,1));
            index_j = kron(ones(4,1),index');
            
            element_K_v1 = reshape(Q1_K_loc',16,1);
    
            K = K + sparse(index_i,index_j,element_K_v1,dim,dim);
                        
            for i = 1 : 4
                rhs(index(i),1) = rhs(index(i),1) + rhs_1_loc(i,1);
            end

        end
        
        P_q1_k= P_q1((k + 1) * (m - 1),:);
        [Q1_K_loc, rhs_1_loc] = q1_local_matrices_rhs(P_q1_k, x_q1, y_q1, n, m, gaussian_nodes, gaussian_weights, number_nodes);
        
        index = zeros(1,2);
        index(1) = P_q1_k(1,1) - ny - 1 - 2 * k;
        index(2) = P_q1_k(1,2) - ny - 3 - 2 * k;
        
        index_boundary = zeros(1,2);
        index_boundary(1) = ny+nx+2+k;
        index_boundary(2) = ny+nx+3+k;
        
        index_i = kron(index',ones(2,1));
        index_j = kron(ones(2,1),index');

        element_K_v1 = reshape(Q1_K_loc([1,2],[1,2])',4,1);
    
        K = K + sparse(index_i,index_j,element_K_v1,dim,dim);
        
        rhs(index(1),1) = rhs(index(1),1) + rhs_1_loc(1,1);
        rhs(index(2),1) = rhs(index(2),1) + rhs_1_loc(2,1);

        index_i_boundary = kron(index',ones(2,1));
        index_j_boundary = kron(ones(2,1),index_boundary');
        
        element_K_v_boundary = reshape(Q1_K_loc([1,2],[4,3])',4,1);

        K_boundary = K_boundary + sparse(index_i_boundary,index_j_boundary,element_K_v_boundary,dim,num_boundary);
        
    end

	% iteration on the right part of the tessellation
    P_q1_k= P_q1((n - 2) * (m - 1) + 1,:);
        
    index = zeros(1,1);
	index(1) = P_q1_k(1,4) - ny - 2 * nx - 1;
    
    index_boundary = zeros(1,3);
	index_boundary(1) = ny+2+nx;
	index_boundary(2) = ny+3+2*nx;
    index_boundary(3) = ny+4+2*nx;

	[Q1_K_loc, rhs_1_loc] = q1_local_matrices_rhs(P_q1_k, x_q1, y_q1, n, m, gaussian_nodes, gaussian_weights, number_nodes);
    
    index_i = index;
    index_j = index;
    
    element_K_v1 = reshape(Q1_K_loc(4,4),1,1);

    K = K + sparse(index_i,index_j,element_K_v1,dim,dim);
    
	rhs(index(1),1) = rhs(index(1),1) + rhs_1_loc(2,1);
    
    index_i_boundary = kron(index',ones(3,1));
    index_j_boundary = kron(ones(1,1),index_boundary');
        
    element_K_v_boundary = reshape(Q1_K_loc(4,1:3)',3,1);

    K_boundary = K_boundary + sparse(index_i_boundary,index_j_boundary,element_K_v_boundary,dim,num_boundary);
    
	for k = 2 : m - 2
        P_q1_k= P_q1((n - 2) * (m - 1) + k,:);
    
    	[Q1_K_loc, rhs_1_loc] = q1_local_matrices_rhs(P_q1_k, x_q1, y_q1, n, m, gaussian_nodes, gaussian_weights, number_nodes);
        
        index = zeros(1,2);
        index(1) = P_q1_k(1,1) - ny - 2 * nx - 1;
        index(2) = P_q1_k(1,4) - ny - 2 * nx - 1;
    
        index_boundary = zeros(1,2);
        index_boundary(1) = ny+2*nx+k+2;
        index_boundary(2) = ny+2*nx+k+3;
        
        index_i = kron(index',ones(2,1));
        index_j = kron(ones(2,1),index');
        
        element_K_v1 = reshape(Q1_K_loc([1,4],[1,4])',4,1);
    
        K = K + sparse(index_i,index_j,element_K_v1,dim,dim);
                
        for i = 1 : 2
            rhs(index(i),1) = rhs(index(i),1) + rhs_1_loc(i,1);
        end
        
        index_i_boundary = kron(index',ones(2,1));
        index_j_boundary = kron(ones(2,1),index_boundary');
        
        element_K_v_boundary = reshape(Q1_K_loc([1,4],[2,3])',4,1);

        K_boundary = K_boundary + sparse(index_i_boundary,index_j_boundary,element_K_v_boundary,dim,num_boundary);

	end
    
    P_q1_k= P_q1((n - 1) * (m - 1),:);
        
    index = zeros(1,1);
	index(1) = P_q1_k(1,1) - ny - 2 * nx - 1;
    
    index_boundary = zeros(1,3);
    index_boundary(1) = ny+2+2*nx;
    index_boundary(2) = 2*(ny+nx+1)+1;
    index_boundary(3) = 2*(ny+nx+2);
    
	[Q1_K_loc, rhs_1_loc] = q1_local_matrices_rhs(P_q1_k, x_q1, y_q1, n, m, gaussian_nodes, gaussian_weights, number_nodes);
    
    index_i = index;
    index_j = index;

    element_K_v1 = reshape(Q1_K_loc(1,1),1,1);

    K = K + sparse(index_i,index_j,element_K_v1,dim,dim);
        
	rhs(index(1),1) = rhs(index(1),1) + rhs_1_loc(1,1);
    
    index_i_boundary = kron(index',ones(3,1));
    index_j_boundary = kron(ones(1,1),index_boundary');
        
    element_K_v_boundary = reshape(Q1_K_loc(1,[4,2,3])',3,1);

    K_boundary = K_boundary + sparse(index_i_boundary,index_j_boundary,element_K_v_boundary,dim,num_boundary);

end