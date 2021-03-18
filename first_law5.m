% without affinity terms
function [F] = first_law5(pi,pf,Ti,Tf,rvi,rvf,rli,rlf,rii,rif)

  global R_d R_v c_pd c_v c_l c_i
    
  Tav = .5*(Ti+Tf);
  rvav = .5*(rvi+rvf);
  rlav = .5*(rli+rlf);
  riav = .5*(rii+rif);
  c = c_pd + rvav*c_v + rlav*c_l + riav*c_i;
    
  %F = (c_pd +rvav*c_v)*(log(Tf) - log(Ti)) - ...
  %    (R_d +rvav*R_v)*(log(pf) - log(pi)) + ...
  %    l_lv(Tav)/Tav*(rvf - rvi);
  F = log(Ti) + 1/c*((R_d + rvav*R_v)*log(pf/pi));

end


