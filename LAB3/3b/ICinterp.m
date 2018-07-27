function [ Uold ] = ICinterp( C, in_cond, totnodes, nodes, elements)

Uold = zeros(totnodes, 1);

for (i = 1:totnodes)
    Uold(i) = in_cond( nodes(i,1), nodes(i,2) );
end

end