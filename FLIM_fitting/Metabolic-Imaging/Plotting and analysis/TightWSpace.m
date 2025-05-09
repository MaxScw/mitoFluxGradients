function TightWSpace(h,Bnd,BndX)
% Simple function to make axis limits almost tight, but leave a spacer
% around the boundary.
% INPUTS:
% -h: axis handle (often 'gca')
% -Bnd: boundary spacer, as fraction of axis dimension (e.g. 0.07)
% -BndX: optional independent bound for X, otherwise, assume equal to Bnd
if ~exist('Bnd') Bnd = 10; end
if ~exist('BndX') BndX = Bnd; end
Ch = get(h,'Children');
X = []; Y = [];
for i = 1:length(Ch)
    if strcmp(Ch(i).Type,'line')|strcmp(Ch(i).Type,'bar')|strcmp(Ch(i).Type,'errorbar')
        X = [X Ch(i).XData];
        Y = [Y Ch(i).YData];
    end
end
Xlms = [min(X) max(X)]; Xrng = Xlms(2)-Xlms(1);
if Xlms(1)==Xlms(2) Xlms(2)=Xlms(1)+.01; end
Ylms = [min(Y) max(Y)]; Yrng = Ylms(2)-Ylms(1);
if Ylms(1)==Ylms(2) Ylms(2)=Ylms(1)+.01; end
set(h,'Ylim',[Ylms(1)-Yrng*Bnd Ylms(2)+Yrng*Bnd]);
set(h,'Xlim',[Xlms(1)-Xrng*BndX Xlms(2)+Xrng*BndX]);
