function [ K ] = Boundary_Dirch_K(K, nelemx, nelemy)

%identify the nodes that are going to be deleted 
nnodesx = nelemx+1;
nnodesy = nelemy+1;
totnodes = (nelemx+1)*(nelemy+1);
del_nodes_2 = [1:1:nelemx+1];
del_nodes_3 = [2*(nelemx+1) : nelemx+1 : totnodes];
del_nodes_1 = [nnodesx+1:nnodesx:totnodes-nnodesx+1];
del_nodes_4 = [totnodes-nnodesx+1:1:totnodes];

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

for (i=del_nodes_1)
    K(i,:) = zeros(1, totnodes);
    K(:,i) = zeros(totnodes,1);
    K(i,i) = 1;
end

for (i=del_nodes_4)
    K(i,:) = zeros(1, totnodes);
    K(:,i) = zeros(totnodes,1);
    K(i,i) = 1;
end

%centered source
k = floor(totnodes/2)+1;
K(k,:) = zeros(1, totnodes);
K(:,k) = zeros(totnodes,1);
K(k,k) = 1;



end

