function h=EBplotyy(x,y)
%INPUT:
%x-independent variable (often Time)
%y-dependent variable
%OUTPUT:
%h-handle of the errorbar graphics object
%
%Save this code in an m-file named EBplotyy.m

s=nanstd(y);
h=errorbar(x,nanmean(y),nanstd(y)./sqrt(size(y,1)),'o','MarkerSize',10);  %NOTE: You might want to change nanmedian to nanmean, nanstd to standard error, etc depending on what you want to plot