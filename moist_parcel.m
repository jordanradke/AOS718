% physical parameters
cd '/Users/jordan/Documents/MATLAB/aos718'

global R_d R_v c_pd c_v c_l c_i eps l_lvt p_ref
R_d     = 287.05;       % gas constant (dry air)
R_v     = 461.5;        % gas constant (vapor)
c_pd    = 1005;         % heat capacity of dry air (pressure cnst)
c_v     = 1870;         % heat capacity of vapor   (pressure cnst)
c_l     = 4190;         % heat capacity of liquid  (pressure cnst)
c_i     = 2106;         % heat capacity of ice
eps     = 1.61;         % R_v/R_d
l_lvt   = 2.5008e6;     % latent heat (vapor -> liquid) at triple point
p_ref   = 1000;         % reference pressure at surface (hectopascals)


% initial conditions 
i = 1;
p(i)   = 1000;   % hPa
T(i)   = 300;    % K
r_v(i) = 15e-3;  % 15 g/kg
r_l(i) = 0;
r_i(i) = 0;

theta(i)   = T(i)*(p_ref/p(i))^(R_d/c_pd);
theta_m(i) = T(i)*(p_ref/p(i))^(R_m(r_v(i))/c_pm(r_v(i),0,0));   
Hl(i) = min(2,H_l(p(i),T(i),r_v(i)));
es(i) = e_star(T(i));
ev(i) = eps^(-1)*r_v(i)*p(i)/(1+r_v(i)*eps);

% pressure increment + final pressure
dp = 1;      % hPa
p1 = 150;    % hPa

% main loop: updates through first law for mixed systems
while p(i) > p1

% compute diagnostic variables
theta(i)   = T(i)*(p_ref/p(i))^(R_d/c_pd);
theta_m(i) = T(i)*(p_ref/p(i))^(R_m(r_v(i))/c_pm(r_v(i),0,0));   
Hl(i) = min(2,H_l(p(i),T(i),r_v(i)));     
es(i) = e_star(T(i));
ev(i) = eps^(-1)*r_v(i)*p(i)/(1+r_v(i)*eps);


% pressure + vapor update straightforwardly
p(i+1)   = p(i)-dp;
r_v(i+1) = r_v(i);

% compute temperature at next step from first law for mixed phase system
% (152) (adapted to this simpler case)
%fun    = @(Tf) first_law1(p(i),p(i+1),T(i),Tf,r_v(i),r_v(i+1));
%T(i+1) = fzero(fun,T(i));  % l.h.s = 0
T(i+1) = exp(first_law1(p(i),p(i+1),T(i),T(i),r_v(i),r_v(i+1)));

i = i+1;
end

theta(i)   = T(i)*(p_ref/p(i))^(R_d/c_pd);
theta_m(i) = T(i)*(p_ref/p(i))^(R_m(r_v(i))/c_pm(r_v(i),0,0));   
Hl(i) = min(2,H_l(p(i),T(i),r_v(i)));    
es(i) = e_star(T(i));
ev(i) = eps^(-1)*r_v(i)*p(i)/(1+r_v(i)*eps);

temp1    = T;
theta1   = theta;
theta_m1 = theta_m;


%% plots
% Temp (Celcius) vs. pressure
%hold on
%figure
%plot(T-273,p)
%set(gca, 'YDir','reverse')
%hold off

% Temp, theta, theta_m, and Hl vs. pressure
hold on
T_ax = -100:1:40;
p_ax = 100:1000;
line(T-273,      p,'Color','r')
line(theta_m-273,p,'Color','b')
line(theta-273,  p,'Color','g')
xlabel('temperature','Interpreter','latex','FontSize',24)
set(gca, 'YDir','reverse')
ax1        = gca;
ax1.XColor = 'r';
ax1.YColor = 'r';

ax1_pos = ax1.Position;
ax2     = axes('Position',ax1_pos,...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none');
line(Hl, p,'Parent',ax2,'Color','k')
line(r_v,p,'Parent',ax2,'Color','c')
set(gca,'YDir','reverse')
title('pressure vs. $T$, $\theta$ and $\theta_m$', 'Interpreter', 'latex', 'FontSize',24);
xlabel('$H_l$','Interpreter','latex','FontSize',24)
ylabel('pressure', 'Interpreter','latex','FontSize',24)
hold off




