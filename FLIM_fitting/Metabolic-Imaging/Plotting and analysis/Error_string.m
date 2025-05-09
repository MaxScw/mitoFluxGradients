function errstr = Error_string(meas,PMerr)
% Simple little function for taking the measurement and uncertainty
% (plus-minus error bars), and outputing a 'plus-minus' string with same
% precision as the leading significant figure of the uncertainty.
% NOTE: getting the 'plus-minus' character is really hard to save to excell
% from Matlab. Instead, I used '_', so just open in Excell and replace-all.
% INPUTS:
% -meas: measurement
% -PMerr: plus-minus error value

% Examples:
% meas = 1.0343;
% PMerr = 0.0029191;
% PMerr = .0014;

[meas PMerr];
% Convention is that if the first digit of error is 1.[something less than 5],
% you include 2 sig figs, otherwise only 1. Or comment out the last
% condition to have 2 sig figs for anything below 2.
if roundsd(PMerr,1,'floor')==getprecision(roundsd(PMerr,1,'floor'))%&roundsd(PMerr,1,'floor')==roundsd(PMerr,1,'round')
    Err1dig = roundsd(PMerr,2);
else
    Err1dig = roundsd(PMerr,1);
end
Prec = getprecision(Err1dig); Prec2 = -log10(Prec);
errstr = [num2str(round(meas,Prec2)) '_' num2str(Err1dig)];

end

function p = getprecision(x)
f = 14-9*isa(x,'single'); % double/single = 14/5 decimal places.
s = sprintf('%.*e',f,x);
v = [f+2:-1:3,1];
s(v) = '0'+diff([0,cumsum(s(v)~='0')>0]);
p = str2double(s);
end