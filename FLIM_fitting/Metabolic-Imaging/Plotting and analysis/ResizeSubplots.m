function ResizeSubplots(h,Bnds,Inter)
% Input handle to a figure and this will resize all the axes to have
% specified spacings. Specify spacings in normalized units (fraction of
% figure
% INPUTS: 
% -h: handle to figure containing subplots
% -Bnds: Boundary spacing. 4 elements - [left, lower, right, upper]
% -Inter: Inter-axis spacing - [x, y]
% h = gcf;
% Bnds = [.05 .05 .03 .04];
% Inter = [.05 .07];


ch = get(h,'Children');
% Clear non-axes children
nonax = [];
for i = 1:length(ch)
    if ~strcmp(ch(i).Type,'axes') nonax = [nonax i]; end
end
ch(nonax) = [];

% Get axis positions
Pos = nan(length(ch),4);
% Sometimes same positions differ by machine precision. Make sure that
% doesn't happen by rounding to the 5th decimal place.
for i = 1:length(ch)
    Pos(i,:) = ch(i).Position;
    ch(i).XLabel.Units = 'normalized';
    ch(i).YLabel.Units = 'normalized';
end 

% Get numbers of rows and columns, sort subplots by rows and columns
% Note, can't really use unique because it catches floating point errors
% that make the same values look different. Have to do manually
Urws = uniqueprc(Pos(:,2),5); Ucls = uniqueprc(Pos(:,1),5);
Nrws = length(Urws); Ncls = length(Ucls); 
for i = 1:Nrws
    for j = 1:Ncls
        sps(i,j) = ch((abs(Pos(:,2)-Urws(i))<10^-5) & (abs(Pos(:,1)-Ucls(j))<10^-5));
    end
end
% Calculate space taken by labels.
YlabW = abs(ch(1).YLabel.Position(1)*Pos(1,3)); % Ylabel width
XlabH = abs(ch(1).XLabel.Position(2)*Pos(1,4)); % Xlabel Height

% Adjust positions 
% XlabW = min(Urws)-(1-(max(Urws)+Pos(1,4))); 
% YlabW = min(Ucls)-(1-(max(Ucls)+Pos(1,3))); 
SpW = (1-Ncls*YlabW-(Ncls-1)*Inter(1)-Bnds(1)-Bnds(3))/Ncls; % subplot axes width
SpH = (1-Nrws*XlabH-(Nrws-1)*Inter(2)-Bnds(2)-Bnds(4))/Nrws; % subplot axes height
Xpos = nan(length(ch),2);
Ypos = nan(length(ch),2);
for i = 1:Nrws
    for j = 1:Ncls
        sps(i,j).Position(1) = Bnds(1)+YlabW+(j-1)*(SpW+YlabW+Inter(1));
        sps(i,j).Position(2) = Bnds(2)+XlabH+(i-1)*(SpH+XlabH+Inter(2));
        sps(i,j).Position(3) = SpW;
        sps(i,j).Position(4) = SpH;
        1;
    end
end 



% Xpos(i,:) = ch(i).XLabel.Position(1:2);
% Ypos(i,:) = ch(i).YLabel.Position(1:2);
    