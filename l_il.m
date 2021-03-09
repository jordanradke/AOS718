% compute latent heat given temperature (Kirchoff's law)
function [L_il] = l_il(T)

    global l_lvt c_i c_l

    L_il = l_lvt + (c_l - c_i)*(T-273.16);

end