function [ Fe ] = Boundary_Dirch_F(F, K, nelemx, nelemy, nodes, gd_2, gd_3, t)
%this function modifies the load vector adding the Dirichlet boudnary
%conditions
 
totnodes = (nelemx+1)*(nelemy+1);
totelem = nelemx+nelemy;
Fe = zeros(totnodes,1);
Fe = F;
%identify the nodes that are going to be deleted 
del_nodes_2 = [1:1:nelemx+1];
del_nodes_3 = [2*(nelemx+1) : nelemx+1 : totnodes];

for (i=del_nodes_2)
    GD_2 = gd_2( nodes(i,1), nodes(i,2), t); 
    Fe = Fe - K(:,i)*GD_2 ;
    
end

for (i=del_nodes_3)
    GD_3 = gd_3( nodes(i,1), nodes(i,2), t);
    Fe = Fe - (K(:, i))*GD_3 ;
    
end

% Now I can substitue the right values 
for (i = del_nodes_2)
    Fe(i) = gd_2( nodes(i,1), nodes(i,2), t);
end
for (i = del_nodes_3)
    Fe(i) = gd_3( nodes(i,1), nodes(i,2), t);
end


end

