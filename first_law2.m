% first law allowing for condensation/evaporation
function [F] = first_law2(pi,pf,Ti,Tf,rvi,rvf)

  global R_d R_v c_pd c_v
    
  Tav = .5*(Ti+Tf);
  rvav = .5*(rvi+rvf);
  c = c_pd + rvav*c_v;
  %F = (c_pd +rvav*c_v)*(log(Tf) - log(Ti)) - ...
  %    (R_d +rvav*R_v)*(log(pf) - log(pi)) + ...
  %    l_lv(Tav)/Tav*(rvf - rvi);
  F = log(Ti) + 1/c*((R_d + rvav*R_v)*log(pf/pi) - l_lv(Tav)/Tav*(rvf-rvi));

end


