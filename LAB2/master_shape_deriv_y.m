function [ fxx ] = master_shape_deriv_y( i, x, y )
% this function returns the value of the i-th basis function derivative wrt to n on
% the master element evaluated in the point (x,y)

if (i==1)
    fx = @(e,n) -(1/4)*(1-e);
    fxx =fx(x,y);
end

if (i==2)
    fx = @(e,n) -(1/4)*(1+e);
    fxx =fx(x,y);
end

if (i==3)
    fx = @(e,n) (1/4)*(1+e);
    fxx =fx(x,y);
end

if (i==4)
    fx = @(e,n) (1/4)*(1-e);
    fxx =fx(x,y);
end




end

