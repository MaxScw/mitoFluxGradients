function y=linear_model(data, p, ~)
    % param = [slope, intercept]
    s = p(1, :);
    intercept = p(2, :);
    y = s.*data + intercept;
end