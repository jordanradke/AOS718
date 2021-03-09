% compute relative humidity given pressure, temp, vapor mixing ratio
function [H] = H_l(p,T,r_v)
    
    H = e_v(p,r_v)/e_star(T);
    
end

% partial pressure w.r.t vapor
function [ev] = e_v(p,r_v)

    global eps
    
    ev = eps^(-1)*r_v*p/(1+r_v*eps);
    
end