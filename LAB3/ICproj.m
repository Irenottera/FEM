function [ Uold ] = ICproj( C, in_cond, totnodes, nodes, elements, t0, nelemx, nelemy)

a = size(elements);
totelem = a(2);
b = zeros (totnodes,1);

for (i = 1:totelem)
    be = LocalstiffF ( i, nodes, elements, in_cond , t0 );
    b = b + be;
end

A = C;
Uold = A\b;

end

