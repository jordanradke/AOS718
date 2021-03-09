% a simplification of the three-phase open system first law when rv is
% constant and no ice or water is involved.
function [F] = first_law1(pi,pf,Ti,Tf,r_v)

    global R_d c_pd c_l
              
    F = c_pd*log(Tf) - R_d*log(p_d(pf,r_v)) + r_v*l_lv(Tf)/Tf - ...
        (c_pd*log(Ti) - R_d*log(p_d(pi,r_v)) + r_v*l_lv(Ti)/Ti) + ...
        r_v*c_l*log(Tf) - r_v*c_l*log(Ti);  % r.h.s. of first law
    
end


