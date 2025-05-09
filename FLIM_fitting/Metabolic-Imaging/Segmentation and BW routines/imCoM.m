function [x,y] = imCoM(mask);
% Simple function for quickly calculating the center of mass of any binary
% image (ie, mask). 

[y,x] = find(mask);
x = mean(x);
y = mean(y);