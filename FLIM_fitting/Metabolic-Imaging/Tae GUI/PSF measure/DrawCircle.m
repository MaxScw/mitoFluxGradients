function h = DrawCircle(xc,yc,radius)
% Making a circle the size of the feature around each feature.
% h = DrawCircle(xc,yc,radius)

if length(radius) == 1
   radius = ones(length(xc),1)*radius; 
end
theta = 0:0.1:2*pi;
h = zeros(length(xc),1);
for i = 1:length(xc)
    x = xc(i) + radius(i)*cos(theta)*2;
    y = yc(i) + radius(i)*sin(theta)*2;
    h(i) = plot(x,y,'-g','LineWidth',1.5);
end