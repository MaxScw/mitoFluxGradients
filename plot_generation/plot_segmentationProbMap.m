%% designation of oxygen ranges to analyse
clear all;
clc;

oxygen_ranges={'run_oxy_l0_h0p1', 'run_oxy_l0p1_h0p2', 'run_oxy_l0p2_h0p3', 'run_oxy_l0p3_h0p4',...
               'run_oxy_l0p4_h0p6', 'run_oxy_l0p6_h0p8',...
               'run_oxy_l0p8_h1', 'run_oxy_l1_h1p2',...
               'run_oxy_l1p2_h1p4', 'run_oxy_l1p4_h1p6',...
               'run_oxy_l1p6_h1p8',...
               'run_oxy_l1p8_h2', 'run_oxy_l2_h2p8',...
               'run_oxy_l2p8_h3p6', 'run_oxy_l3p6_h4p4',...
               'run_oxy_l4p4_h5p2'};
num_oxygen_ranges={[0, 0.1], ...
                   [0.1, 0.2], [0.2, 0.3], [0.3, 0.4],...
                   [0.4, 0.6], [0.6, 0.8], [0.8, 1.0], [1.0, 1.2],...
                   [1.2, 1.4], [1.4, 1.6], [1.6, 1.8], [1.8, 2.0],...
                   [2.0, 2.8], [2.8, 3.6], [3.6, 4.4],... 
                   [4.4, 5.2]};

load('../data/exp_solubilities.mat')
load('../data/exp_temperatures.mat')

temp_resolved = false;
pp_path = '../data/EXP_published/postprocessing_results/';
plot_path = '../data/EXP_published/plots/';
dpath = "../data/EXP_published/FLIM_fitting_results/flim_structs";

% temp_resolved = true;
% pp_path = '/home/mx/mitoFluxGradients/data/EXP_temperatureResolved/postprocessing_results/';
% plot_path = '/home/mx/mitoFluxGradients/data/EXP_temperatureResolved/plots';
% temp_ind = 1;
if temp_resolved==true
    temp_string = '_T'+string(temperature(temp_ind))+'C';
    dpath = '../data/EXP_temperatureResolved/FLIM_fitting_results/temp_'+string(temperature(temp_ind))+'C';
else 
    temp_string = '';
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