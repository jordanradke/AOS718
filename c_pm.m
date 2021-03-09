% compute moist heat capacity given vapor mixing ratio
function [cpm] = c_pm(r_v)

    global c_pd c_v
    
    cpm = c_pd + r_v*c_v;
    
end