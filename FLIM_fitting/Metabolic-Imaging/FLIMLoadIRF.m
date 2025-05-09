function IRFstruct = FLIMLoadIRF(IRFpath,adcratio)
% Load IRF, automatically pick off fit, noise, and peak bounds from a irf so that
% you don't have to manually adjust when acquisition settings are changed.
% Finally, clean up IRF
% Inputs:
% -irf: 1D decay vector for IRF
% -adcratio: if you want the bounds for a 256 bin decay, enter 16
% (4096/256). The bounds will be output for the 256 decays.
% Outputs:
% -fst/end: fit start and end
% -nst/end: noise start and end
% -pst/end: irf peak start and end
% -bst/bend: background start and end

% VERSIONS:
% 2015-05-23: Change noise range to be the last 10th of the IRF instead of
% at the beginning

% clear all;
% IRFpath = 'Z:\Tim\2015-04-08,14 Cembs Temp variation\2015-04-14 More Temp Tests\irf_40xobj_ex750nm_no_em_int30sec.sdt';

sdt = bh_readsetup(IRFpath);
dt = sdt.SP_TAC_TC*10^9;
irf = double(bh_getdatablock(sdt,1));
% Try to pick off knee, but if there are not enough photons, it's not
% possible to detect that spike. Thus, just set knee to about 1/5 in
% (usually where the spike comes). Doesn't really matter, because the
% spike position will be determined later from the IRF.
L = length(irf);
IRFstruct.irforig = irf;

% Weird BH glitch. Sometimes irf has '1' for the last element of the irf
% array. Set all '1's to 0
irf(irf==1)=0;
sm=round(L/100);
smdecay=(smooth(irf,sm));
ddecay = diff(smdecay);
% Take stdev of the initial noise. Find the start index
ind1 = find(irf); %That's where you start getting signal
% Get the noise level (stddev) from the last 10th of the decay
rng = (ind1(end)-round(L/10)):(ind1(end));
stdev = std(ddecay(rng));
%Finally, find the index where the derivative starts to shoot up. Set
%thresh for 7*stdev
knee = find(ddecay>7*stdev); knee = knee(1);

% That was the hard part. Now set all the bounds to reasonable values
nonzero = find(irf);
fst = nonzero(1); fend = nonzero(end);
nst = fst; nend = round(knee-(knee-fst)/10);
pst = round(knee-(knee-fst)/10); pend = round(knee + L/17);

if exist('adcratio')
    fst = ceil(fst/adcratio);
    fend = floor(fend/adcratio);
    nst = ceil(nst/adcratio);
    nend = floor(nend/adcratio);
end

% Clean up IRF

% Keep these noise bounds for purposes of fitting the decays. For decays,
% the noise is better calculated right before the intenisity spike.
% HOWEVER, the 750nm IRF w/ no filter has a scattering artifact there, so
% for cleaning up the IRF, calculate the background from the end section of
% the array.
backgr = mean(irf((fend-round((fend-fst)/10)):fend));
irf([nst:pst,pend:end]) = 0;
irf(pst:pend) = max(irf(pst:pend) - backgr,0);

dt = sdt.SP_TAC_TC*10^9;
time = (1:length(irf))'*dt;
save([IRFpath(1:end-4) '.mat'],'irf','time');

% Renormalize to max = 5000 counts, just to make it convenient to plot IRFS
% against each other. The normalization doesn't matter because it gets
% renormalized in lm_decay_model.m to match the data's total photon counts.
irf = irf/max(irf)*5*1E3;

IRFstruct.fst = fst;
IRFstruct.fend = fend;
IRFstruct.nst = nst;
IRFstruct.nend = nend;
IRFstruct.irf = irf;
IRFstruct.time = time;
IRFstruct = orderfields(IRFstruct);
