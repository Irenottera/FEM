function [ Uold ] = ICproj( C, in_cond, totnodes, nodes, elements)


b = zeros (totelem);
for (i = 1:totelem)
    be = LocalstiffF ( e, nodes, elements, q, t0 );
    b = b + be;
end

A = C;

Uold = A\b;



end

