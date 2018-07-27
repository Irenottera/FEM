function [ C ] = Cmatrix( nelemx, nelemy, nodes, elements)
%this fucntion creates the C matrix - time invariant 

totnodes = (nelemx+1)*(nelemy+1);
totelem = nelemx*nelemy;

C = zeros(totnodes, totnodes);

for e = 1:totelem
    
        Ce = zeros(totnodes, totnodes);
        cur_nodes = nodes(elements(e, :)',:);
        
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
        
        ii = 1;
        jj = 1;
        
        for i = elements(e, :)
            jj = 1;
            for j = elements(e, :)
                
                %this is TO CHECK!! 
                k1 = PHI(1,ii)*PHI(1,jj);
                k2 = PHI(2,ii)*PHI(2,jj);
                k3 = PHI(3,ii)*PHI(3,jj);
                k4 = PHI(4,ii)*PHI(4,jj);

                Ce(i,j) = w1*k1+w2*k2+w3*k3+w4*k4 ;
                
                jj = jj+1;
            end
            ii = ii+1;
        end
        
     C = C + Ce;   
        
end






end

