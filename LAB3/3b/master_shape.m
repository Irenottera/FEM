function [ fxx ] = master_shape( i, x, y )
% this function returns the value of the i-th basis function on the master element in (x,y)
% i is the function index i=1:4
% x, y are the point in which the function needs to be evaluated 

if (i==1)
    fx = @(e,n) (e-1)*(n-1)*(1/4) ;
    fxx =fx(x,y);
end

if (i==2)
    fx = @(e,n) (e+1)*(1-n)*(1/4) ;
    fxx =fx(x,y);
end

if (i==3)
    fx = @(e,n) (e+1)*(n+1)*(1/4) ;
    fxx =fx(x,y);
end

if (i==4)
    fx = @(e,n) (1-e)*(n+1)*(1/4) ;
    fxx =fx(x,y);
end



end

