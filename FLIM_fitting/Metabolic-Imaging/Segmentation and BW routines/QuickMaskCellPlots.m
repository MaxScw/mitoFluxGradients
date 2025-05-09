%%
figure
m1 = zeros(512,512);
for i =1:length(masks) 
    m1 = m1 + masks{i};
end
imshow(m1)
for i =1:length(masks)
    [Cy,Cx] = imCoM(masks{i});
    text(Cy,Cx,num2str(i),'fontsize',30,'color','r')
end


close all
m1 = zeros(512,512);
for i =1:length(LastFrMasks)
    m1 = m1 + LastFrMasks{i};
end
imshow(m1)
for i =1:length(LastFrMasks)
    [Cy,Cx] = imCoM(LastFrMasks{i});
    text(Cy,Cx,num2str(i),'fontsize',30,'color','r')
end
