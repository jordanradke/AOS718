% compute dry partial pressure from total pressure + vapor
function [pd] = p_d(p,r_v)

    global eps
    
    pd = p*(1+ r_v/eps)^(-1);

end