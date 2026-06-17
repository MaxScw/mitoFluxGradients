%% load setup file for file paths and data
clc
clear all
temp_resolved = false;
run("setup.m")

%% define path for FLIM structs
dpath = data_path+"EXP_published/FLIM_fitting_results/flim_structs";
if temp_resolved==true
    dpath = '../data/EXP_temperatureResolved/FLIM_fitting_results/temp_'+string(temperature(temp_ind))+'C';
end

%% read-in of data

load_index = 7;
filename_cellular_kn=dir(dpath + '/*');

load(dpath + '/' + filename_cellular_kn(load_index).name);

mitoprob = flim_struct(2).mitoprob;
cytoprob = flim_struct(2).cytoprob;
backprob = flim_struct(2).backprob;

nadhint = flim_struct.img;

%%
% thresholding like in original paper 
thresh_mitoprob = mitoprob;
thresh_mitoprob(thresh_mitoprob<=0.7) = 0;
thresh_mitoprob(thresh_mitoprob>0.7) = 1;

thresh_cytoprob = cytoprob;
thresh_cytoprob(thresh_cytoprob<=0.7) = 0;
thresh_cytoprob(thresh_cytoprob>0.7) = 1;

thresh_backprob = backprob;
thresh_backprob(thresh_backprob<=0.7) = 0;
thresh_backprob(thresh_backprob>0.7) = 1;

xStart = 230%1;
xEnd = 430%220;
yStart = 65%250;
yEnd = 263%470;

fig = figure('Renderer', 'painters', 'Position', [10 10 500 500]);
imshow(nadhint)
hold on
clim([0, 15])
xlim([xStart xEnd])
ylim([yStart yEnd])

PixelSize = 0.42;
Scalebar_length = 20;
x_location = 380; 
y_location = 245;
quiver(x_location,y_location,Scalebar_length/PixelSize,0, ...
       'ShowArrowHead','off', Color=[1 1 1], LineWidth=4)
text(x_location+7, y_location+8, "20$\mu$m", ...
    Color=[1 1 1], Interpreter="latex",FontSize=20)

C = round([656 332]/2) ;   % center of circle 
th = linspace(0,2*pi) ; 

mito_color = uisetcolor([1 0 1]);

% draw concentric circle 
for r = linspace(5, 80, 10)
    x = C(1)+r*cos(th) ; 
    y = C(2)+r*sin(th) ; 
    plot(x,y, 'Color',mito_color, LineStyle='--', LineWidth=2)
end

savefig(string(plot_path)+'NADH_intMap'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'NADH_intMap'+string(temp_string)+'.png')
%saveas(fig, string(plot_path)+'NADH_intMap'+string(temp_string)+'.svg')
export_fig(fig, string(plot_path)+'NADH_intMap'+string(temp_string)+'.eps')
%%
fig2 = figure('Renderer', 'painters', 'Position', [10 10 500 500]);
imshow(thresh_mitoprob - thresh_backprob)
%imshow(thresh_mitoprob - thresh_cytoprob)
clim([-1 1])
mymap = [
    0 0 0
    1 1 1
    1 0 1
    ];
colormap(mymap)

xlim([xStart xEnd])
ylim([yStart yEnd])
savefig(string(plot_path)+'NADH_intMap_segment'+string(temp_string)+'.fig')
saveas(fig2, string(plot_path)+'NADH_intMap_segment'+string(temp_string)+'.png')
export_fig(fig2, string(plot_path)+'NADH_intMap_segment'+string(temp_string)+'.eps')