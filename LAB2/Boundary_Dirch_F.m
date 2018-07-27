function [ Fe ] = Boundary_Dirch_F(F, K, nelemx, nelemy, nodes, elements, gd_2, gd_4)
%this function modifies the load vector adding the Dirichlet boudnary
%conditions
 
totnodes = (nelemx+1)*(nelemy+1);
Fe = zeros(totnodes,1);
Fe = F;
%identify the nodes that are going to be deleted 
del_nodes_2=[1:1:nelemx+1];
del_nodes_4=[(totnodes-nelemy):1:totnodes];

%vector of values of the g_d functions on the nodes that are being deleted
GD_2 = 0;
GD_4 = 0;

for (i=del_nodes_2)
    GD_2 = gd_2( nodes(i,1), nodes(i,2)); 
    Fe = Fe - K(:,i)*GD_2 ;
    
end

for (i=del_nodes_4)
    GD_4 = gd_4( nodes(i,1), nodes(i,2));
    Fe = Fe - (K(:, i))*GD_4 ;
    
end

end

