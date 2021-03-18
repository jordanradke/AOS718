% compute latent heat given temperature (Kirchoff's law)
% l_lvt = latent heat at triple point
function [L_lv] = l_lv(T)

    global l_lvt c_v c_l

    L_lv = l_lvt + (c_v - c_l)*(T-273.16);

end
