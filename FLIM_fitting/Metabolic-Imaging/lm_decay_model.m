function y = lm_decay_model(time,param,consts,irf)

% % TEST - save vars for inputs for fast testing and checking
% clear all; load('C:\Users\Tim\Documents\Academic - Research\Data\temp_synthIRF.mat')
% clear all; load('C:\Users\Tim\Documents\Academic - Research\Data\temp_wIRF.mat')

%% Versions:
% 2019-05-22: updated to allow for generating a synthetic IRF, if no IRF is
% provided to the function. Also implemented a faster fft convolution.
% NOTE: I went ahead and implemented it, but it didn't seem to give very
% close result to the physical IRF for the test data. Do a more careful
% comparison when it's useful, but for now, just make sure to use measured
% IRFs, as usual. 

% normalized to the number of photons
%

% MODEL:
% param(1) : shift of decay model from IRF (usually ranges from -10 to 10)
% Functions
% 1-expo: P = A*exp(-t/tau1) + (1-A)
% P = param(2)*exp(-t/param(3))+(1-param(2))
% param(2): fractional amp of expo decay(s). Fraction of photons coming
%           from fluorescence signal vs dark noise background. (1-A) = bg
% param(3): lifetime of the first decay
%
% 2-expo: P = A*[(1-f2)*exp(-t/tau1)+ f2*exp(-t/tau2)] + (1-A)
% P = param(2)*( (1-param(4))*exp(-t/param(3))+param(4)*exp(-t/(param(5))) )+(1-param(2));
% param(4) : fraction of the second decay
% param(5) : lifetime of second decay
%
% 3-expo: P = A(1-f2-f3)*exp(-t/tau1)+ A*f2*exp(-t/tau2) + A*f3*exp(-t/tau3) + (1-A)
% P = param(2)*( (1-param(4)-param(6))*exp(-t/param(3))+param(4)*exp(-t/(param(5)))+param(6)*exp(-t/(param(7))) )+(1-param(2));
% param(6) : fraction of the third decay
% param(7) : lifetime of third decay
% 
% consts = [nexpo,counts(seg),fit_start,fit_end,dt_irf,LaserT,GaussM_init]
% nexpo: number of exponentials in your model
% counts: Total number of photons in the region of interest
% dt_irf: time bin size for IRF
% LaserT: Period of laser pulses.

nexpo = consts(1);
counts = consts(2);
fit_start = consts(3);
fit_end = consts(4);
shift = round(param(1));
% NOTE: If an IRF is provided, 'shift' is the amount the IRF is shifted.
% For synthetic IRF, 'shift' is the total offset of the mean from zero.
len_ratio = 16;  %oversampling rate

% Define time axis (small bin size, like IRF)
dt_irf = (time(2)-time(1))/len_ratio;
T = consts(6); % Laser rep rate defines full time window
t = (0:dt_irf:T)';
    
% IRF: if time_irf and irf are not provided, or they are not arrays, 
% irf is assumed to be a gaussian with 
% mean = param(1)  and std = param(end);
% I.e. the gaussian std is a fit parameter, added to the end of the array

% INFER IRF?
if ~exist('irf')|irf==-1
    %length check
    if length(param) ~= (nexpo+1)*2
        error('length(param) should be 2*(nexpo+1)');
    end
    
    % Define IRF mean and std, and convert from bins to time units
    GaussMeanInitShift = consts(7)/length(time)*T;
    irf_mean = param(1)/length(t)*T + GaussMeanInitShift;
    irf_std = param(end)/length(t)*T;
    
    % Inferred IRF
    irf = exp(-(t-irf_mean).^2/(2*irf_std^2));
    irf = irf/sum(irf);
else
    irf = double(irf/(sum(irf)));
    
    %length check
    if length(param) ~= 2*nexpo+1
        error('length(param) should be 2*nexpo+1');
    end
    shift = round(param(1));
    irf = mycircshift(irf,shift);
end

if nexpo == 1
    decay_model = param(2)*exp(-t/param(3))+(1-param(2));
elseif nexpo == 2
    decay_model = param(2)*( (1-param(4))*exp(-t/param(3))+param(4)*exp(-t/(param(5))) )+(1-param(2));
elseif nexpo == 3
    decay_model = param(2)*( (1-param(4)-param(6))*exp(-t/param(3))+param(4)*exp(-t/(param(5)))+param(6)*exp(-t/(param(7))) )+(1-param(2));
end

decay_model(end+1:end+length(t)) = decay_model(1:length(t));

% Fourier Transform of model
nfft = 2^nextpow2(2*length(t)+length(irf)-1); % Precisions of fourier trans
Fmodel = fft(decay_model,nfft).*fft(irf,nfft);
% 
Nt = length(t);
convmodel = real(ifft(Fmodel));
% Alt - Matlab
% convmodel = conv(irf,decay_model)*dt_irf;
convmodel = convmodel(Nt+1:Nt+length(time)*len_ratio);

new_y = reshape(convmodel,len_ratio,length(time));
new_y = sum(new_y,1);
y = new_y';
% Note, ends of model aren't truncated to fit window, but in lm_matx in
% lm.m, only that section is compared to the data.

%normalize to the number of total counts
y = y/sum(y(fit_start:fit_end))*counts;

% % test with (comment/uncomment): 
% % Matlab convolution, same result and a little slower
% yml = conv(irf,decay_model)*dt_irf;
% yml = yml(Nt+1:Nt+length(time)*len_ratio);
% yml = reshape(yml,round(len_ratio),length(time));
% yml = sum(yml,1);
% yml = yml/sum(yml(fit_start:fit_end))*counts;
% subplot(1,4,1); plot(decay_model); title('decay model'); 
% subplot(1,4,2); plot(irf); title('irf'); 
% subplot(1,4,3); plot(y); title('irf convnfft w model'); 
% subplot(1,4,4); plot(yml); title('irf conv(matlab) w model'); 
% set(gcf,'position',[100 100 1400 400])
