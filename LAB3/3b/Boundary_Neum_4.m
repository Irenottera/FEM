function [ Ge ] = Boundary_Neum_4(e, nelemx, nelemy, nodes,nnodesx, nnodesy, elements, gn_4,  t)

totelem = nelemx*nelemy;
totnodes = (nelemx+1)*(nelemy+1);

% I find the elements that ly on the right boundary 
    %cur_elem = zeros(nelemy, 4);
    cur_nodes = zeros (4, 2);
    cur_nodes = nodes(elements(e, :)',:);
      
    %2 integration points  - remember it's a 1-d integral
    int_points = zeros (2,2);
    int_points(1, :) = [ -1/sqrt(3) , 1];
    int_points(2, :) = [1/sqrt(3) , 1];
    %weight functions 
    w1 = 1;
    w2 = 1;

    %basis functions on the 2 integration points 
    PHI=zeros(2,4);
    X = zeros(2,1);
    Y = zeros(2,1);
    GN = zeros(2,1);
    DXDN = zeros(2,1);
    DYDN = zeros(2,1);
    J = zeros (2,1);
    Ge = zeros(totnodes, 1);
    for (i = 1:2)
        for (j = 1:4)
            PHI(i,j) = master_shape(j, int_points(i, 1), int_points(i,2));
            X(i) = X(i) + cur_nodes(j,1)*master_shape(j,int_points(i,1),int_points(i,2));
            Y(i) = Y(i) + cur_nodes(j,2)*master_shape(j,int_points(i,1),int_points(i,2));
            DXDN(i) = DXDN(i) + cur_nodes(j,1)*master_shape_deriv_y(j,int_points(i,1),int_points(i,2)); 
            DYDN(i) = DYDN(i) + cur_nodes(j,2)*master_shape_deriv_y(j,int_points(i,1),int_points(i,2));
        end
    end

    for (i=1:2)
        GN(i) = gn_4(X(i), Y(i), t);
    end

    for (i=1:2)
        J(i) = sqrt( ( DXDN(i) )^2 + ( DYDN(i) )^2 );
    end

    ii=1;
    for (i=elements(e, :))
        Ge(i) = GN(1)*PHI(1,ii)*J(1)*w1 + GN(2)*PHI(2,ii)*J(2)*w2;
        ii = ii+1;
    end
    


end

