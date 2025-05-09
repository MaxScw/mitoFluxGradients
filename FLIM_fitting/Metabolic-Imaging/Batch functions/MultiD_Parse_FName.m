function [T,ch,Z] = MultiD_Parse_FName(fname)
% Short function to look at a file name of form
% 'sdt_000000000_NADH_002.sdt' and extract the uMan macro T, Z, and channel
% information
dashes = strfind(fname,'_');
T = str2num(fname(dashes(end-2)+1:dashes(end-1)-1));
ch = fname(dashes(end-1)+1:dashes(end)-1);
Z = str2num(fname(dashes(end)+1:dashes(end)+3));