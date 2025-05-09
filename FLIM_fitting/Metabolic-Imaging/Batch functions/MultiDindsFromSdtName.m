function [Tind,Chind,Zind] = MultiDindsFromSdtName(fname)

dashes = strfind(fname,'_');
Tind = str2num(fname(dashes(1)+1:dashes(2)-1))+1;
chan = fname(dashes(2)+1:dashes(3)-1);
if strcmp(chan,'NADH')
    Chind = 1;
elseif strcmp(chan,'FAD')
    Chind = 2;
else
    Chind = 3;
end
Zind = str2num(fname(dashes(3)+1:end-4))+1;

