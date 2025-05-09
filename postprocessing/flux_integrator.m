function concentration=flux_integrator(data, param, c)
    % param = [c_star, D]
    c_star = param(1);
    D = param(2);
    jox = c;

    function dy_dr=rhs_first_order_system(r, y)
        c = y(1);
        chi = y(2);
        dc_dr = chi;
        dchi_dr = interp1(data, jox, r)./D - 2.*chi./r;
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
