function CircMask = Make_Circ_Mask(xdim,ydim,Cx,Cy,r)
% Quick function to make a circular mask of radius, r, centered around
% [Cx,Cy]

if isempty('Cx')|Cx==-1 Cx = round(xdim/2); end
if isempty('Cy')|Cy==-1 Cy = round(ydim/2); end

CircMask = zeros(ydim,xdim);
[Xm,Ym] = meshgrid((1:xdim)-Cx,(1:ydim)-Cy);
CircMask(sqrt(Xm.^2+Ym.^2)<r)=1;