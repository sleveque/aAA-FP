function [gaussian_nodes, gaussian_weights] = gaussian_nodes_and_weights(number_nodes)

	if number_nodes == 2
        gaussian_nodes = zeros(1,2);
        gaussian_weights = zeros(1,2);
    
        gaussian_nodes(1,1) = -sqrt(3)/3;
        gaussian_nodes(1,2) = sqrt(3)/3;
    
        gaussian_weights(1,1) = 1;
        gaussian_weights(1,2) = 1;
        
	elseif number_nodes == 3
        gaussian_nodes = zeros(1,3);
        gaussian_weights = zeros(1,3);
    
        gaussian_nodes(1,1) = -sqrt(0.6);
        gaussian_nodes(1,2) = 0;
        gaussian_nodes(1,3) = sqrt(0.6);
    
        gaussian_weights(1,1) = 5/9;
        gaussian_weights(1,2) = 8/9;
        gaussian_weights(1,3) = 5/9;

	end

end
