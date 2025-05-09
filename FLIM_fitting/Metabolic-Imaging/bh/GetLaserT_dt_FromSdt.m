function [dt,laserT] = GetLaserT_dt_FromSdt(sdtsetup)
% Calculate time bin size of the irf, and get laser rep rate (period, actually)
% Used in lm_decay_model.m.

% dt time bin size
range = sdtsetup.SP_TAC_R*10^9;
gain = double(sdtsetup.SP_TAC_G);
resol = double(sdtsetup.SP_ADC_RE);
dt = range/(gain*resol); 

%get laser rep period and period
%default value (in ns)
stdmeas = bh_getmeasdesc(sdtsetup,1);
laserT = 1/stdmeas.min_sync_rate*10^9;

