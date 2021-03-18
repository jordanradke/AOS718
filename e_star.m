% formula for saturation vapor pressure over liquid
function [estar] = e_star(T)

    % convert to celcius
    Tc = T - 273.15;
    
    a = 17.2693882;
    b = 35.86;
    estar = 6.1078*exp(a*(T-273.16)/(T-b));        % thermo handbook
    
    
    %if Tc >= 0
    %    estar = .61078*exp(17.27*Tc/(T+237.3));    % Wikipedia
    %else
    %    estar = .61078*exp(21.875*Tc/(T+265.5));
    %end
    
    %estar = .611*10^(7.5*Tc/(Tc +237.3));          % AMS
    
    % From Huang 2018
    %if Tc > 0
    %    estar = exp(34.494-4924.99/(Tc+237.1))/(Tc+105)^(1.57);
    %else
    %    estar = exp(43.494 - 6545.8/(Tc+868)^2);
    %end
    
end