function concentration=theory_flux_integrator(data, param, const, param_profiles)
    % param = [c_star, D]
    c_star = param(1);
    D = param(2);
    if param_profiles==false
        vmax_expParams = const(1:3);
        km_linParams = const(4:5);
        corr_fact = const(6:end, :);
    else
        vmax_profile = const(1, :);
        km_profile = const(2, :);
    end

    function dy_dr=rhs_first_order_system(r, y)
        c = y(1);
        chi = y(2);
        dc_dr = chi;
        if param_profiles==false
        dchi_dr = ((interp1(data, corr_fact, r).*exp_model(r, vmax_expParams).*c)./(linear_model(r, km_linParams) + c))./D - 2.*chi./r;
        else
        dchi_dr = ((interp1(data, vmax_profile, r).*c)./(interp1(data, km_profile, r) + c))./D - 2.*chi./r;
        end
        dy_dr = [dc_dr
                 dchi_dr];
    end

    function res=boundary_condition(y0, yR)
        c0 = y0(1);
        chi0 = y0(2);
        cR = yR(1);
        dy_dr0 = rhs_first_order_system(0, [c0, chi0]);
        dc_dr0 = dy_dr0(1);
        res = [dc_dr0-1e-3
               cR-c_star];
    end

    function y = guess(r)
        y = [c_star
             1e-6];
    end

    r_range = data;
    c_solinit = bvpinit(r_range, @guess);
    c_sol = bvp5c(@rhs_first_order_system, @boundary_condition, c_solinit);
    c_data = deval(c_sol, r_range);
    
    concentration = c_data;
end
