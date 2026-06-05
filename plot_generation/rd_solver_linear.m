function flux=rd_solver_linear(data, param, c)
    % param = [R, A, lambda]
    R = param(1);
    A = param(2);
    lambda = param(3);
    flux = A.*(R./data).*(sinh(data./lambda))./(sinh(R/lambda));
end