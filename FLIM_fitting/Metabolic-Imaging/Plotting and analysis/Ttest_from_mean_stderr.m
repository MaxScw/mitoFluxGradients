function p = Ttest_from_mean_stderr(mns,stderrs,ns)
% Enter the sample means, standard errors, and sample sized for two
% samples, and this function calculates the t-test p-values for you. 
% Note: written for 2018 paper with Denny, in order to compare embryo
% metabolic parameter values at different time points, taking stderr from
% the error on the measurement (e.g. uncertainty in NADH fraction engaged,
% returned by lm fitting routine).
% Inputs:
%  -mns: 2-element array with sample means from both samples
%  -stderrs: standard errors for both samples
%  -ns: sample number for both samples (for FLIM params, #of photons)

% NOTE: OBS. I did z-tests instead for that paper. See:
% 'C:\Users\Tim\Documents\Academic - Research\Publication Materials\2018 Denny Intro to MetIm\Figures\Fig4\RevisedAnalysis\FLIM_mdpars_B2_B1_merge_Timeplots_phists1CB_3Zs_fin.m'
% z-tests = z=(mu1-mu2)/sqrt(se1^2+se2^2), then (1-normcdf(abs(z1c(rw,p,ch))))*2;

Amean =mns(1); 
Asd=stderrs(1);
Bmean=mns(2); 
Bsd=stderrs(2);
v = ns(1)+ns(2)-2;
tval = (Amean-Bmean) / sqrt((Asd^2+Bsd^2));       % Calculate T-Statistic
% tval = (Amean-Bmean) / sqrt((Asd^2/ns(1)+Bsd^2/ns(2)));       % If we were going to use std dev instead of std err.
tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution
tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;              % 1-tailed t-distribution
tprob = 1-[tdist2T(tval,v)  tdist1T(tval,v)]