function [fst, fend, nst, nend, pst, pend] = FLIMGetFitBoundsFromIRF(irf,adcratio)
% OBSOLETE: Simplify IRF loading by creating a new function called
% 'FLIMLoadIRF', and calculate the bounds within that.
% Automatically pick off fit, noise, and peak bounds from a irf so that
% you don't have to manually adjust when acquisition settings are changed.
% Inputs:
% -irf: 1D decay vector for IRF
% -adcratio: if you want the bounds for a 256 bin decay, enter 16
% (4096/256). The bounds will be output for the 256 decays.
% Outputs:
% -fst/end: fit start and end
% -nst/end: noise start and end
% -pst/end: irf peak start and end
% -bst/bend: background start and end

% load('C:\Users\Tim\Desktop\test\dec.mat')


% Try to pick off knee, but if there are not enough photons, it's not
% possible to detect that spike. Thus, just set knee to about 1/5 in
% (usually where the spike comes). Doesn't really matter, because the
% spike position will be determined later from the IRF.
L = length(irf);

% Weird BH glitch. Sometimes irf has '1' for the last element of the irf
% array. Set all '1's to 0
irf(irf==1)=0;
sm=round(L/100);
smdecay=(smooth(irf,sm));
ddecay = diff(smdecay);
% Take stdev of the initial noise. Find the start index
ind1 = find(irf); ind1 = ind1(1); %That's where you start getting signal
% Define the range as the 80% of noise between on-point and signal spike
rng = (ind1+round(L/20)):(ind1+round(2*L/20));
stdev = std(ddecay(rng));
%Finally, find the index where the derivative starts to shoot up. Set
%thresh for 7*stdev
knee = find(ddecay(rng(end):end)>7*stdev)+rng(end); knee = knee(1);

% That was the hard part. Now set all the bounds to reasonable values
nonzero = find(irf);
fst = nonzero(1); fend = nonzero(end);
nst = fst; nend = round(knee-(knee-fst)/10);
pst = round(knee-(knee-fst)/10); pend = round(knee + L/7);

if exist('adcratio')
    fst = ceil(fst/adcratio);
    fend = floor(fend/adcratio);
    nst = ceil(nst/adcratio);
    nend = floor(nend/adcratio);
end