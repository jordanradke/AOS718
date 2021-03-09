% Teten's formula for saturation vapor pressure over ice
function [estar] = e_starice(T)

    a = 21.8745584;
    b = 7.66;
    estar = 6.1078*exp(a*(T-273.16)/(T-b));
    
end