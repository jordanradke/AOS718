% compute moist gas constant given vapor mixing ratio
function [Rm] = R_m(r_v)

    global R_d R_v
    
    Rm = R_d + r_v*R_v;
    
end