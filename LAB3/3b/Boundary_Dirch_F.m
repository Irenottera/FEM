function [ Fe ] = Boundary_Dirch_F(F, K, nelemx, nelemy, nodes, gd_1, gd_2, gd_3, gd_4, t)
%this function modifies the load vector adding the Dirichlet boudnary
%conditions

nnodesx=nelemx+1;
totnodes = (nelemx+1)*(nelemy+1);
totelem = nelemx+nelemy;
Fe = zeros(totnodes,1);
Fe = F;

%identify the nodes that are going to be deleted 
del_nodes_2 = [1:1:nelemx+1];
del_nodes_3 = [2*(nelemx+1) : nelemx+1 : totnodes];
del_nodes_1 = [nnodesx+1:nnodesx:totnodes-nnodesx+1];
del_nodes_4 = [totnodes-nnodesx+1:1:totnodes];

%vector of values of the g_d functions on the nodes that are being deleted

for (i=del_nodes_1)
    GD_1 = gd_1( nodes(i,1), nodes(i,2), t); 
    Fe = Fe - K(:,i)*GD_1 ;
    
end

for (i=del_nodes_2)
    GD_2 = gd_2( nodes(i,1), nodes(i,2), t); 
    Fe = Fe - K(:,i)*GD_2 ;
    
end

for (i=del_nodes_3)
    GD_3 = gd_3( nodes(i,1), nodes(i,2), t);
    Fe = Fe - (K(:, i))*GD_3 ;
    
end
for (i=del_nodes_4)
    GD_4 = gd_4( nodes(i,1), nodes(i,2), t); 
    Fe = Fe - K(:,i)*GD_4 ;
    
end

% Now I can substitue the right values 
for (i = del_nodes_1)
    Fe(i) = gd_1( nodes(i,1), nodes(i,2), t);
end
for (i = del_nodes_2)
    Fe(i) = gd_2( nodes(i,1), nodes(i,2), t);
end
for (i = del_nodes_3)
    Fe(i) = gd_3( nodes(i,1), nodes(i,2), t);
end
for (i = del_nodes_4)
    Fe(i) = gd_4( nodes(i,1), nodes(i,2), t);
end

%centered source 
k = floor(totnodes/2)+1;
if (t<=0.01)
    GD = 1+100*t; 
    Fe = Fe - K(:,k)*GD ;
    Fe(k) = 1+100*t;
elseif(t>0.01)
    GD = 2; 
    Fe = Fe - K(:,k)*GD ;
    Fe(k) = 2;
end


end

