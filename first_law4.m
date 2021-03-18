% first law allowing for condensation/evaporation// condensate removed
% immediately
function [F] = first_law4(pi,pf,Ti,Tf,rvi,rvf,rli,rlf,rii,rif)

  global R_d R_v c_pd c_v
    
  Tav = .5*(Ti+Tf);
  rvav = .5*(rvi+rvf);
  rlav = .5*(rli+rlf);
  riav = .5*(rii+rif);
  c = c_pd + rvav*c_v;
  
  %F = (c_pd +rvav*c_v)*(log(Tf) - log(Ti)) - ...
  %    (R_d +rvav*R_v)*(log(pf) - log(pi)) + ...
  %    l_lv(Tav)/Tav*(rvf - rvi);
  F = log(Ti) + 1/c*((R_d + rvav*R_v)*log(pf/pi) - ...
      l_lv(Tav)/Tav*(rvf-rvi)) + l_il(Tav)/Tav*(rif-rii);

end




