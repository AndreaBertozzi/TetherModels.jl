"""
    get_analytic_catenary(filename)

Loads the analytic catenary curve for the 2D catenary example

# Arguments
- filename: the filename of the mat file to read

# Returns
- x_cat: x coordinates of the catenary curve
- y_cat: x coordinates of the catenary curve
"""
function get_analytic_catenary(filename)
    vars        = matread(filename)
    vars        = get(vars, "analytic_catenary", 0)
    x_cat       = vec(get(vars, "x", 0))
    y_cat       = vec(get(vars, "y", 0))
    return x_cat, y_cat
end

function transformFromOtoW(windDirection_rad,vec_O)
    M_WO = [cos(windDirection_rad) sin(windDirection_rad) 0; sin(windDirection_rad) -cos(windDirection_rad) 0;  0 0 -1]
    vec_W = M_WO*vec_O
    return vec_W
end

function transformFromWtoO(windDirection_rad,vec_W)
    M_OW = [cos(windDirection_rad) sin(windDirection_rad) 0; sin(windDirection_rad) -cos(windDirection_rad) 0; 0 0 -1]
    vec_O = M_OW*vec_W
    return vec_O
end

function calculate_rho_at_height(h, set)
    if h <=0
        rho_at_height = set.rho_0
    else
        rho_at_height = set.rho_0*exp(-h/set.h_p)
    end
    return rho_at_height
end

end
