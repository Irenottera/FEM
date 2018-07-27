function [ elements ] = RELEM( totelem, nelemx )
%This function creates a matrix where in every i-th row we have the number
%of the 4 nodes that describe the i-th element of the mesh. 

elements = zeros(totelem, 4);
for (i = 1: totelem) %number of the element
        if (mod(i, nelemx)~=0)
            elements(i, 1) =  i+ floor(i/nelemx);
            elements(i, 2) =  elements(i, 1) + 1 ;
            elements(i, 3) =  i+ floor(i/nelemx) + nelemx +2;
            elements(i, 4) =  i+ floor(i/nelemx) + nelemx +1;
        elseif ( mod(i, nelemx)==0)
            elements(i, 1) =  i+ floor(i/nelemx)-1;
            elements(i, 2) =  elements(i, 1)+1 ;
            elements(i, 3) =  i+ floor(i/nelemx) + nelemx +1;
            elements(i, 4) =  i+ floor(i/nelemx) + nelemx ;
        end
    
end

end

