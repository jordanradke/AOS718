% Teten's formula for saturation vapor pressure over liquid
function [estar] = e_star(T)

    a = 17.2693882;
    b = 35.86;
    estar = 6.1078*exp(a*(T-273.16)/(T-b));
    
end