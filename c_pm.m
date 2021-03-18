% compute moist heat capacity given vapor mixing ratio
function [cpm] = c_pm(r_v,r_l,r_i)

    global c_pd c_v c_l c_i
    
    cpm = c_pd + r_v*c_v + r_l*c_l + r_i*c_i;
    
end