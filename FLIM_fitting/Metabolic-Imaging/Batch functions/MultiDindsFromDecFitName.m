function [Pstr,Mind] = MultiDindsFromDecFitName(fname)

dashes = strfind(fname,'_');
Pstr = fname(dashes(1)+4:dashes(2)-1);
Mind = str2num(fname(dashes(2)+4:dashes(3)-1));

