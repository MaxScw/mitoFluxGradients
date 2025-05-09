function IJpath = MatIJpathconv(path)
% Simple function for converting paths in Windows/Matlab format into the
% format that ImageJ macros recognize. Simply replaces all '\'s with '\\'s

% path = 'C:\Users\Tim\Documents\Academic - Research\Data\SegTest';

if path(end)~='\' path = [path '\']; end
slashes = strfind(path,'\');
IJpath = path;
for i = length(slashes):-1:1
    if i == length(slashes)
        IJpath = [IJpath '\'];
    else
        IJpath = [IJpath(1:slashes(i)) '\' IJpath(slashes(i)+1:end)];
    end
end