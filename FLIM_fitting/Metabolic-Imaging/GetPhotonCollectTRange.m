function AcqRng = GetPhotonCollectTRange(sdtfile)
% Uses metadata in an sdt file and extracts the exact time range during
% which photons were being collected. Units are absolute, datenum time
% (days)

sdt = bh_readsetup(sdtfile);
meas = bh_getmeasdesc(sdt,1);
numscans = double(meas.hist_fida_points);

collectend = datenum([sdt.Date ' ' sdt.Time]);
% stopt = double(meas.stop_time);
cellectduration = double(meas.fcs_end_time);
collectst = collectend - cellectduration/86400;

AcqRng = [collectst collectend];