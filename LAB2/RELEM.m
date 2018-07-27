function [ elements ] = RELEM( totelem, nelemx )

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

