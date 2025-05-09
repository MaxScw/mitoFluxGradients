function MultiD_Position_Plot_prototype(file)
% Give program the saved positions from a uManager acquisition, and this
% plots the positions with numbers so you can double check the positions
% were entered correctly.

% clear all;
% file = 'Z:\Lab\Tim\2016-09-01 Emre collab 2\Dish1_oldyoung.pos';

mat =TextScanToColumn(file);
% mat = VarName1;
Pinds = [];
for i = 1:length(mat)
    if strfind(mat{i},'Pos')
        ind = strfind(mat{i},'Pos');
        qts = strfind(mat{i},'"');
        pos{i} = mat{i}(ind+3:qts(end)-1);
        Pinds = [Pinds i];
        xcoord = mat{i-20};
        x(i) = str2num(xcoord(6:end));
        ycoord = mat{i-13};
        y(i) = str2num(ycoord(6:end));
    end
end
emp = find(x==0);
pos(emp)= [];x(emp)=[];y(emp)=[];
plot(x,y,'.')
text(x,y,pos)