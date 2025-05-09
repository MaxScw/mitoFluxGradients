function [p_srt,sg_srt] = SrtFLIMparamsByTaus(p_fit,sigma_p,corr)
% Reorder fit params in order of increasing lifetimes.
% param(1) : shift of decay model from IRF (usually ranges from -10 to 10)
% Functions
% 1-expo: P = A*exp(-t/tau1) + (1-A)
% P = param(2)*exp(-t/param(3))+(1-param(2))
% param(2): fractional amp of expo decay(s). Fraction of photons coming
%           from fluorescence signal vs dark noise background. (1-A) = bg

% After that, subsequent pairs indicate [f_n,tau_n].
% E.g. p_fit params for a 3 exp should be
% [shf, A, shortest tau, f2(of middle tau), middle tau, f3 (longest tau), logest tau]
% then f1 = 1-f2-f3;
% clear all;
% load('C:\Users\Tim\Documents\Academic - Research\Data\REWRITE2EXP\p_fit_outoforder.mat')
% 2-exp test - swap order
% p_fit(4) = 1-p_fit(4); % Should be frac of longer LT
% p_fit = p_fit([1 2 5 4 3]); sigma_p = sigma_p([1 2 5 4 3]);
% 3-exp test:
% p_fit=[p_fit;.11;.1]; sigma_p = [sigma_p;.123;.124];
% p_fit = p_fit([1 2 3 6 7 4 5]);
% p_fit=[p_fit;.1;2]; sigma_p = [sigma_p;.123;.124];

% Case of synthetic IRF, then p_fit has an extra param on the end for the
% gaussian width. Remove it, sort pars, then re-add
corr = 1; % Dummy val
SynIRFBool = 0;
if mod(length(p_fit),2)==0
    SynIRFBool = 1;
    GaussW = [p_fit(end) sigma_p(end) corr(end)];
    p_fit(end) = []; sigma_p(end) = []; corr(end) = [];
end

if length(p_fit) == 3 % 1-exp
    % [shf, A, only tau]
    p_srt = p_fit;
    sg_srt = sigma_p;
elseif length(p_fit) == 5 % 2-exp
    % [shf, A, shortest tau, f2(of longest tau), longest tau]
    if p_fit(3)>p_fit(5) % if wrong order, switch
        SrtInds = [1 2 5 4 3];
        p_srt = p_fit(SrtInds);
        p_srt(4) = 1-p_srt(4); % Should be frac of longer LT
        sg_srt = sigma_p(SrtInds);
%         for i = 1:length(p_fit)
%             for j=1:length(p_fit)
%                 corr_srt(i,j) = corr(SrtInds(i),SrtInds(j));
%             end
%         end
%         % Change sign of fraction correlations, since fraction changed
%         corr_srt(4,:)=-corr_srt(4,:); corr_srt(:,4)=-corr_srt(:,4);
    else
        p_srt = p_fit;
    sg_srt = sigma_p;
    end
elseif length(p_fit) == 7 % 3-exp
    clear SrtInds
    % [shf, A, shortest tau, f2(of middle tau), middle tau, f3 (of longest tau), logest tau]
    LTinds = 3:2:length(p_fit);
    LTs = p_fit(LTinds); [a,LTord]=sort(LTs);
    % Make a matrix of [fractions, taus, fraction sigma_ps, tau sigma_ps]
    FracTauPairs(1:3,1) = [1-p_fit(4)-p_fit(6);p_fit(4);p_fit(6)];
    FracTauPairs(1:3,2) = p_fit(3:2:7);
    FracTauPairs(1:3,3) = [sqrt(sigma_p(4)^2+sigma_p(6)^2);sigma_p(4);sigma_p(6)];
    FracTauPairs(1:3,4) = sigma_p(3:2:7);
    FracTauPairs(1:3,5) = [1 0 0]; % Show which was the previous t1
    FracTauPairs = FracTauPairs(LTord,:);
    % Make sorted p's. Construct SrtInds as we go
    p_srt(1:2,1) = p_fit(1:2); SrtInds(1:2) = [1 2];
    % Fill in LTs
    p_srt(3:2:length(p_fit)) = FracTauPairs(:,2); 
    SrtInds(3:2:length(p_fit)) = LTinds(LTord);
    % Fill in fractions f2 and f3
    p_srt([4 6]) = FracTauPairs(2:3,1);
    % If fractions have had change of variables, we need to keep track, and
    % list indices of previous fractions that are being combined to form
    % new one. E.g. f3'=f1=1-f2-f3, then keep track of [-4 -6] to
    % calculate new correlation coefficients.
    SrtInds(2,:) = 0; Frinds = [4 6];
    for i = 2:3
        if FracTauPairs(i,5)==1
            SrtInds(1:2,Frinds(i-1)) = [-4; -6];
        else
            SrtInds(1,Frinds(i-1)) = Frinds(i-1);
        end
    end
    
    % Make sorted sigma_ps
    sg_srt(1:2,1) = sigma_p(1:2);
    % Fill in LTs
    sg_srt(3:2:length(p_fit)) = FracTauPairs(:,4);
    % Fill in fractions f2 and f3
    sg_srt([4 6]) = FracTauPairs(2:3,3);
    
%     % Resort timcorrelation matrix and recalculate corr coefs for change of
%     % fraction variables (what a pain)
%     % NOTE: I finished this, but I can't really be sure if it's correct,
%     so don't use if for now. Get someone else to check it.
%     for i = 1:length(p_fit)
%         for j=1:length(p_fit)
%             if i == 4 & j ==4
%                 1;
%             end
%             corr_srt(i,j) = sign(SrtInds(1,i))*sign(SrtInds(1,j))*corr(abs(SrtInds(1,i)),abs(SrtInds(1,j)))+...
%                 sign(SrtInds(2,i))*sign(SrtInds(2,j))*corr(abs(SrtInds(1,i)),abs(SrtInds(1,j)));
%         end
%     end
    
    % Compare
    %     [p_fit sigma_p p_srt sg_srt]
end

if SynIRFBool
    p_srt = [p_srt; GaussW(1)]; sg_srt = [sg_srt; GaussW(2)]; corr = [corr GaussW(3)];
end