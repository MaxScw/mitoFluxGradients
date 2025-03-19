function [log_prob]=norm_eval(par)

mu=par(1);
sig=par(2);
var=par(3);

log_prob=log(1./sqrt(2.*pi.*sig^2).*exp(-(var-mu)^2./(2.*sig^2)));

end

