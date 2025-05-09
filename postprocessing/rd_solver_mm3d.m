function flux=rd_solver_mm3d(data, param, ~)
    % param = [c_star, v_max, k_m, D]
    c_star = param(1);
    v_max = param(2);
    k_m = param(3);
    D = param(4);
    

    %k_m = (v_max/joxR - 1)*c_star;
    %v_max = (k_m/c_star + 1)*joxR;

    function dy_dr=rhs_first_order_system(r, y)
        c = y(1);
        chi = y(2);
        dc_dr = chi;
        dchi_dr = v_max.*c./(D.*(c + k_m)) - 2.*chi./r;
        dy_dr = [dc_dr
                 dchi_dr];
    end

    function res=boundary_condition(y0, yR)
        c0 = y0(1);
        chi0 = y0(2);
        cR = yR(1);
        dy_dr0 = rhs_first_order_system(0, [c0, chi0]);
        dc_dr0 = dy_dr0(1);
        res = [dc_dr0
               cR-c_star];
    end

    function y = guess(r)
        % R = max(r);
        % y = [10.*R.*sinh(r./R)./r
        %      10.*R.*cosh(r./R)./r];
        y = [1
             1];
    end

    r_range = data;
    c_solinit = bvpinit(r_range, @guess);
    c_sol = bvp4c(@rhs_first_order_system, @boundary_condition, c_solinit);
    c_data = deval(c_sol, r_range);
    
    flux = c_data(1, :).*v_max./(c_data(1, :) + k_m); 
    flux = flux';
end
