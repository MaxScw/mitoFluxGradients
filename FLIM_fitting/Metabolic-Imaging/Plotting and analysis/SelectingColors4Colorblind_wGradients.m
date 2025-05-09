
load('colorblind_colormap.mat')
col = colorblind;
close all
for i = 1:size(col,1)
    plot(1:10,(1:10).*i,'color',col(i,:),'LineWidth',2); hold on;
    text(10,10*i,num2str(i));
end
col1 = col(3,:);
col2 = col(6,:);
col3 = col(12,:);
% col1 = [.5 0 0];
% col2 = [1 .7 .2];
% col3 = [0 0 .6];
n1 = 2; n2 = 3;
df1 = (col2-col1)./(n1-1);
df2 = (col3-col2)./(n2-1);
figure;
for i = 1:n1
    gradcols(i,:) = col1+df1.*(i-1);
    plot(1:10,(1:10).*i,'color',col1+df1.*(i-1),'LineWidth',2); hold on;
end
for i = n1+1:n1+n2
    gradcols(i,:) = col2+df2.*(i-n1-1);
    plot(1:10,(1:10).*i,'color',col2+df2.*(i-n1-1),'LineWidth',2); hold on;
end
gradcols(i+1,:) = col3*.7;
plot(1:10,(1:10).*(i+1),'color',col3*.7,'LineWidth',2); hold on;
% save('GradCols.mat','gradcols')