

function [updir subdr] = UpOneDir(base)
% Outputs the directory one level up from the one entered.

if base(end)~='\'
    base = [base '\'];
end

Slashes = find(base=='\');

updir = base(1:Slashes(end-1));
subdr = base(Slashes(end-1)+1:end-1);
