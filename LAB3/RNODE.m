function [ nodes ] = RNODE( nnodesx, nnodesy, x1, y1, xn, yn )
%This function creates a matrix with ascissa and ordiate of every mesh
%node.

nelemx = nnodesx - 1;
nelemy = nnodesy - 1;
%definition of increments = side length of each element 
dx = (xn-x1)/nelemx;
dy = (yn-y1)/nelemy;
totnodes = nnodesx*nnodesy;
%vector of nodes coordinates
nodes = zeros(totnodes, 2);

for (i=1:totnodes)
    
    if (mod(i, nnodesx) ~= 0 ) 
        nodes(i, 1)= (mod(i, nnodesx) -1)*dx ;
        nodes(i, 2)= (floor(i/nnodesx))*dy;
    elseif  (mod(i, nnodesx) == 0)
        nodes(i, 1)= xn ;
        nodes(i, 2)= nodes(i-1, 2);
    end
    
end

end

