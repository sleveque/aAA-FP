function P_q1 = q1_connectivity_matrix(n, m)

    % n: number of intervals on the x-axis => n+1 number of points on the
    % x-axis
    % m: number of intervals on the y-axis => m+1 number of points on the
    % y-axis

    P_q1 = zeros(4, n*m);
    
    for i = 0 : n-1
        for j = 1 : m
            P_q1(1, i*m+j)= i*(m+1)+j;
            P_q1(2, i*m+j)= (i+1)*(m+1)+j;
            P_q1(3, i*m+j)= (i+1)*(m+1)+j+1;
            P_q1(4, i*m+j)= i*(m+1)+j+1;
        end
    end
    
    P_q1 = P_q1';

end