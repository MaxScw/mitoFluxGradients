function y=exp_model(data, p, ~)
    % param = [slope, intercept]
    amp = p(1, :);
    dec = p(2, :);
    offset = p(3, :);
    y = amp.*exp(data.*dec) + offset;
end