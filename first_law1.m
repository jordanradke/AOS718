% a simplification of the three-phase open system first law when rv is
% constant and no ice or water is involved.
% based on (152) in thermo handbook
%function [F] = first_law1(pi,pf,Ti,Tf,r_v)

 %   global R_d c_pd c_l
              
  %  F = c_pd*log(Tf) - R_d*log(p_d(pf,r_v)) + r_v*l_lv(Tf)/Tf - ...
  %      (c_pd*log(Ti) - R_d*log(p_d(pi,r_v)) + r_v*l_lv(Ti)/Ti) + ...
  %      r_v*c_l*log(Tf) - r_v*c_l*log(Ti);  % r.h.s. of first law
    
%end

% first law allowing for condensation/evaporation
% (170)
function [F] = first_law1(pi,pf,Ti,Tf,rvi,rvf)

  global R_d R_v c_pd c_v
    
  Tav = .5*(Ti+Tf);
  rvav = .5*(rvi+rvf);
  c = c_pd + rvav*c_v;
  %F = (c_pd +rvav*c_v)*(log(Tf) - log(Ti)) - ...
  %    (R_d +rvav*R_v)*(log(pf) - log(pi)) + ...
  %    l_lv(Tav)/Tav*(rvf - rvi);
  F = log(Ti) + 1/c*((R_d + rvav*R_v)*log(pf/pi));
      %l_lv(Tav)/Tav*(rvf-rvi) + l_il(Tav)/Tav*(rif - rii));
end




