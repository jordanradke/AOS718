cd '/Users/jordan/Documents/MATLAB/aos718'

% physical parameters
global eps 


% initial conditions 
i = 1;
p(i)   = 1000;   % hPa
T(i)   = 300;    % K
r_v(i) = 15e-3;  % 15 g/kg
r_l(i) = 0;
r_i(i) = 0;
r_t    = r_v(i) + r_l(i) + r_i(i);

theta(i)   = T(i)*(p_ref/p(i))^(R_d/c_pd);
theta_m(i) = T(i)*(p_ref/p(i))^(R_m(r_v(i))/c_pm(r_v(i),0,0));   
Hl(i) = H_l(p(i),T(i),r_v(i));       


% pressure increment + final pressure
dp = 1;      % hPa
p1 = 150;    % hPa

% main loop: updates using first law for mixed systems
while p(i) > p1
    
% first: compute diagnostic variables from last step.
theta(i)   = T(i)*(p_ref/p(i))^(R_d/c_pd);
theta_m(i) = T(i)*(p_ref/p(i))^(R_m(r_v(i))/c_pm(r_v(i),0,0));   
Hl(i) = H_l(p(i),T(i),r_v(i));            

% pressure updates straightforwardly
p(i+1)   = p(i)-dp;

%% iteratively compute temperature and mixing ratios
err  = .01;
j    = 0;
J = 100;

% initial guesses for vapor and liquid mixing ratios
rv_old = r_v(i);
rl_old = r_l(i);
T_old  = T(i);

rv_new = rv_old;
rl_new = rl_old;
T_new  = T_old+1; % so you can enter the loop!

while j < J && abs(T_new-T_old) > err
    % compute new temp given vapor/liquid mixing ratios
    %fun    = @(T) first_law2(p(i),p(i+1),T(i),T,r_v(i),rv_old);
    %T_new  = fzero(fun,T_old);
    T_new = exp(first_law4(p(i),p(i+1),T(i),T_old,r_v(i),rv_old,0,0,0,0));
    
    % compute vapor and liquid mixing ratios given new temp
    r_sat  = eps*e_star(T_new)/(p(i+1) - e_star(T_new));     % saturation ratio
    rv_new = min(r_sat,rv_old);
    
    % enforce conservation of r_t
    %dr_l   = rl_new - rl_old;
    %rv_new = rv_old - dr_l;     % dr_v = -dr_l
    
    % update
    T_old  = T_new;
    rv_old = rv_new;
    rl_old = rl_new;
    j = j+1;
end

%% assign T, p, mixing ratios + update loop
T(i+1) = T_new;
r_v(i+1) = rv_new;
r_l(i+1) = rl_new;

i = i+1;

end

theta(i)   = T(i)*(p_ref/p(i))^(R_d/c_pd);
theta_m(i) = T(i)*(p_ref/p(i))^(R_m(r_v(i))/c_pm(r_v(i),0,0));   
Hl(i) = H_l(p(i),T(i),r_v(i));      

temp4    = T;
theta4   = theta;
theta_m4 = theta_m;


%% plots 
% Temp (Celcius) vs. pressure
hold on
figure
plot(T-273,p)
set(gca, 'YDir','reverse')
hold off

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
%line(r_v/r_t,p,'Parent',ax2,'Color','c')
set(gca,'YDir','reverse')
title('pressure vs. $T$, $\theta$ and $\theta_m$', 'Interpreter', 'latex', 'FontSize',24);
xlabel('$H_l$','Interpreter','latex','FontSize',24)
ylabel('pressure', 'Interpreter','latex','FontSize',24)
hold off


%% comparison of problems 1-5

% temperature in C
hold on
plot(temp1-273, p, 'r')
plot(temp2-273, p, 'b')
plot(temp3-273, p, 'g')
plot(temp4-273, p, '--r')
plot(temp5-273, p, '--b')
set(gca, 'YDir','reverse')
title('pressure vs. temperature', 'Interpreter', 'latex', 'FontSize',24);
xlabel('temperature','Interpreter','latex','FontSize',24)
ylabel('pressure', 'Interpreter','latex','FontSize',24)
hold off

% theta vs pressure
hold on
plot(theta1, p, 'r')
plot(theta2, p, 'b')
plot(theta3, p, 'g')
plot(theta4, p, '--r')
plot(theta5, p, '--b')
set(gca, 'YDir','reverse')
title('pressure vs. $\theta$', 'Interpreter', 'latex', 'FontSize',24);
xlabel('$\theta$','Interpreter','latex','FontSize',24)
ylabel('pressure', 'Interpreter','latex','FontSize',24)
hold off

% theta_m vs pressure
hold on
plot(theta_m1, p, 'r')
plot(theta_m2, p, 'b')
plot(theta_m3, p, 'g')
plot(theta_m4, p, '--r')
plot(theta_m5, p, '--b')
set(gca, 'YDir','reverse')
title('pressure vs. $\theta_m$', 'Interpreter', 'latex', 'FontSize',24);
xlabel('$\theta_m$','Interpreter','latex','FontSize',24)
ylabel('pressure', 'Interpreter','latex','FontSize',24)
hold off



