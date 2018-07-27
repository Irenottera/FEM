function [ Uold ] = IC( C, in_cond, totnodes, nodes, elements)

M = C;
b = zeros(totnodes, 1);

for e=1:totelem
    
    cur_nodes = nodes(elements(e, :)',:);
    
    for j = cur_nodes
        %integration points 
        int_points = zeros (4 , 2);
        int_points(1,1)= -1/sqrt(3);
        int_points(1,2)= -1/sqrt(3);
        int_points(2,1)= 1/sqrt(3);
        int_points(2,2)= -1/sqrt(3);
        int_points(3,1)= 1/sqrt(3);
        int_points(3,2)= 1/sqrt(3);
        int_points(4,1)= -1/sqrt(3);
        int_points(4,2)= 1/sqrt(3);

        %weight functions
        w1 = 1;
        w2 = 1;
        w3 = 1;
        w4 = 1;
        
        %basis functions at integration points 
        PHI = zeros(4,4);
        DXDE = zeros(4,1);
        DXDN = zeros(4,1);
        DYDE = zeros(4,1);
        DYDN = zeros(4,1);
        J = zeros(4,1);
        
        for i=1:4
            for j=1:4
            PHI(i,j) = master_shape(j, int_points(i,1), int_points(i,2));
            DXDE(i) = DXDE(i) + cur_nodes(j,1)*master_shape_deriv_x(j,int_points(i,1),int_points(i,2));
            DXDN(i) = DXDN(i) + cur_nodes(j,1)*master_shape_deriv_y(j,int_points(i,1),int_points(i,2)); 
            DYDE(i) = DYDE(i) + cur_nodes(j,2)*master_shape_deriv_x(j,int_points(i,1),int_points(i,2));
            DYDN(i) = DYDN(i) + cur_nodes(j,2)*master_shape_deriv_y(j,int_points(i,1),int_points(i,2)); 
            end
        end
        
        for i=1:4
            J(i) = DXDE(i)*DYDN(i)-DXDN(i)*DYDE(i);
        end
        
        w1 = w1*J(1);
        w2 = w2*J(2);
        w3 = w3*J(3);
        w4 = w4*J(4);
        
        k1 = u1
        k2 = 
        k3 = 
        k4 = 
        
        b(j) = w1*k1+w2*k2+w3*k3+w4*k4;
        
    end
    
end

Uold = M\b;


end

