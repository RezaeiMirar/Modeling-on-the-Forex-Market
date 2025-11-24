function [a1,a0] = regression(x,y)
n = length(x);
sumx = sum(x);
sumy = sum(y);
sumxx = sum(x.*x);
sumxy = sum(x.*y);
den = n*sumxx - sumx^2;
a1 = (n*sumxy - sumx*sumy)/den;
a0 = (sumxx*sumy - sumxy*sumx)/den;
l = zeros(n,1);
for i=1:n,
    l(i) = a1*x(i) + a0;
end
   

