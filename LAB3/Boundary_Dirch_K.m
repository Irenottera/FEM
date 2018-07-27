function [ K ] = Boundary_Dirch_K(K, nelemx, nelemy, nodes, gd_2, gd_3, t)

%identify the nodes that are going to be deleted 
nnodesx = nelemx+1;
nnodesy = nelemy+1;
totnodes = (nelemx+1)*(nelemy+1);
del_nodes_2 = [1:1:nelemx+1];
del_nodes_3 = [2*(nelemx+1) : nelemx+1 : totnodes];

for (i=del_nodes_2)
    K(i,:) = zeros(1, totnodes);
    K(:,i) = zeros(totnodes,1);
    K(i,i) = 1;
end

for (i=del_nodes_3)
    K(i,:) = zeros(1, totnodes);
    K(:,i) = zeros(totnodes,1);
    K(i,i) = 1;
end


end

