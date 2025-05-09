function un = uniqueprc(arr,pr)
% 'unique' function sometimes sees two equal numbers as different due to
% floating point errors. This considers numbers the same if they are equal
% to below a certain decimal precition.
% INPUTS:
% -arr: input array
% -pr: precision in decimal places

% clear all;
% load('C:\Users\Tim\Downloads\un.mat')
% arr = Pos(:,2); pr = 5;

un = unique(arr);

Rem = [];
for i = 1:length(un)
    for j = 1:length(un)
        if j>i
            if abs(un(i)-un(j))<10^(-pr)
                EssEq(i,j) = 1;
                Rem = [Rem j];
            end
        end
    end
end
Rem = unique(Rem);
un(Rem)=[];

