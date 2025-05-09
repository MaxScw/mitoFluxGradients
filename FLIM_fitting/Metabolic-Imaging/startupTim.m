% Simple script for setting all the figure defaults to display more nicely

if length(find(get(0,'defaultfigurecolor')==1))<3 % First check whether startup has already been run
    
    get(0,'Factory');
    set(0,'defaultfigurecolor',[1 1 1]);
    set(0,'defaultAxesFontSize',12);
    set(0,'DefaultAxesFontSize',13);
    set(0,'DefaultLegendFontSize',12);
    set(0,'DefaultFigurePaperPositionMode','auto')
    
    format short g;
    
end