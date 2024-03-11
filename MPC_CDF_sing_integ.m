clc;
clear;
close all;
addpath("C:\Users\Sajad\Documents\Casadi")
%% Setup and Parameters
x0 = [10; 10];%Initial Condition
xf = [0; 0]; % X final
time_total = 10;%Time for the total steps, equal to tsim
dt = 0.01;
Q = 10*diag([10,10]);
P_weight = 100*diag([10,10]);
R = 1*diag([1]);
N = 20;


%Variables range as below
xmin = [-inf; -inf];
xmax = -xmin;%[10; 10; 10; 100; 100; 100; 100; 100];

umin = [-inf];
umax = -umin;

%First Obstacle
x_obs_1 = 7;
y_obs_1 = 7.3;
r_obs_1 = 1;

%Second Obstacle
x_obs_2 = 4.1;
y_obs_2 = 4.1;
r_obs_2 = 1.5;
%%
import casadi.*
% Symbols/expressions
x = SX.sym('x');
y = SX.sym('y');

states = [x; y];
n_states = length(states);

%Control Inputs
u = SX.sym('u');
controls = [u];
n_controls = length(controls);

% Right hand side
 %---- Parameters

 A = [0 0;
      0 0];
 B = [1;1];

 rhs = A*states + B*controls

% nonlinear mapping function f(x,u)
f = Function('f',{states,controls},{rhs}); 
% Decision variables (controls)
U = SX.sym('U',n_controls,N); 
% parameters (which include at the initial state of the robot and the reference state)
P = SX.sym('P',n_states + n_states);
% A vector that represents the states over the optimization problem.
X = SX.sym('X',n_states,(N+1));

obj = 0; % Objective function
g = [];  % constraints vector

st  = X(:,1); % initial state
g = [g;st-P(1:n_states)]; % initial condition constraints
for k = 1:N
    st = X(:,k);
    con = U(:,k);
    obj = obj+(st-P(n_states+1:2*n_states))'*Q*(st-P(n_states+1:2*n_states)) + con'*R*con; % calculate obj
    st_next = X(:,k+1);
    f_value = f(st,con);
    st_next_euler = st+ (dt*f_value);
    g = [g;st_next-st_next_euler]; % compute constraints
end



%Add Terminal Cost
k = N+1;
st = X(:,k);
obj = obj+(st-P(n_states+1:2*n_states))'*P_weight*(st-P(n_states+1:2*n_states)); % calculate obj

%Adding density function
% tau = SX.sym('tau');
% y = if_else(tau > 0, exp(-1/tau), 0);
% f = Function('f', {tau}, {y});
% f_bar = f(tau)/(f(tau)+f(1-tau));
% Sk = (r_obs+2)^2;
% P_lyap = 0.5*eye(4);
% for k = 1:N
%     st = X(:,k);
%     con = U(:,k);
%     hk = (st(1)-x_obs_1)^2 + (st(2)-y_obs_1)^2 + (st(3)-z_obs_1)^2 - r_obs_1^2;
%     Sk = (st(1)-x_obs_1)^2 + (st(2)-y_obs_1)^2 + (st(3)-z_obs_1)^2 - (r_obs_1+1)^2 ;
%     temp1 = hk / (hk-Sk);
%     f_bar1 = if_else(temp1 > 0, exp(-1/temp1), 0)/(if_else(temp1 > 0, exp(-1/temp1), 0) + if_else(1-temp1 > 0, exp(-1/(1-temp1)), 0));
%     Phi = f_bar1;
%     V = ([st(1) st(2) st(3) st(4)] *P_lyap*[st(1); st(2); st(3); st(4)]);
%     rho = Phi/V;
% 
%     st_next = X(:,k+1);
%     % con_next = U(:,k+1);
%     hk_next = (st_next(1)-x_obs_1)^2 + (st_next(2)-y_obs_1)^2 + (st_next(3)-z_obs_1)^2 - r_obs_1^2;
%     Sk_next = (st_next(1)-x_obs_1)^2 + (st_next(2)-y_obs_1)^2 + (st_next(3)-z_obs_1)^2 - (r_obs_1+1)^2 ;
%     temp2 = hk_next/(hk_next-Sk_next);
%     f_bar2 = if_else(temp2 > 0, exp(-1/temp2), 0)/(if_else(temp2 > 0, exp(-1/temp2), 0) + if_else(1-temp2 > 0, exp(-1/(1-temp2)), 0));
%     Phi_next = f_bar2;
%     V_next = [st_next(1) st_next(2) st_next(3) st_next(4)] *P_lyap*[st_next(1); st_next(2); st_next(3); st_next(4)];
%     rho_next = Phi_next/V_next;
% 
%     f_value = f1(st);
%     f_value_next = f1(st_next);
%     C1 = (f_value_next(1)-f_value(1))/(st_next(1)-st(1)) + (f_value_next(2)-f_value(2))/(st_next(2)-st(2))...
%         + (f_value_next(3)-f_value(3))/(st_next(3)-st(3)) + (f_value_next(4)-f_value(4))/(st_next(4)-st(4))...
%         + (f_value_next(5)-f_value(5))/(st_next(5)-st(5)) + (f_value_next(6)-f_value(6))/(st_next(6)-st(6))...
%         + (f_value_next(7)-f_value(7))/(st_next(7)-st(7)) + (f_value_next(8)-f_value(8))/(st_next(8)-st(8));
% 
%     C1 = C1*rho;
% 
%     g_value = g1(st);
%     g_value_next = g1(st_next);
%     C2 = 0;
%     for i=1:4
%         temp = (g_value_next(1,i) - g_value(1,i)/(st(1)-st_next(1))) + (g_value_next(2,i) - g_value(2,i)/(st(2)-st_next(2)))...
%               + (g_value_next(3,i) - g_value(3,i)/(st(3)-st_next(3))) + (g_value_next(4,i) - g_value(4,i)/(st(4)-st_next(4)))...
%               + (g_value_next(5,i) - g_value(5,i)/(st(5)-st_next(5))) + (g_value_next(6,i) - g_value(6,i)/(st(6)-st_next(6)))...
%               + (g_value_next(7,i) - g_value(7,i)/(st(7)-st_next(7))) + (g_value_next(8,i) - g_value(8,i)/(st(8)-st_next(8)));
%         temp = temp * rho * con(i);
%         C2 = C2 + temp;
%     end
% 
%     g = [g; rho_next - rho + C2 + C1];
%     % g = [g; rho_next ];
% 
% end

% 
% P_lyap = 0.5*eye(4);
% for k = 1:N
%     st = X(:,k);
%     hk = (st(1)-x_obs_2)^2 + (st(2)-y_obs_2)^2 + (st(3)-z_obs_2)^2 - r_obs_2^2;
%     Sk = (st(1)-x_obs_2)^2 + (st(2)-y_obs_2)^2 + (st(3)-z_obs_2)^2 - (r_obs_2+1)^2 ;
%     temp1 = hk / (hk-Sk);
%     f_bar1 = if_else(temp1 > 0, exp(-1/temp1), 0)/(if_else(temp1 > 0, exp(-1/temp1), 0) + if_else(1-temp1 > 0, exp(-1/(1-temp1)), 0));
%     Phi = f_bar1;
%     V = ([st(1) st(2) st(3) st(4)] *P_lyap*[st(1); st(2); st(3); st(4)]);
%     rho = Phi/V;
% 
%     st_next = X(:,k+1);
%     hk_next = (st_next(1)-x_obs_2)^2 + (st_next(2)-y_obs_2)^2 + (st_next(3)-z_obs_2)^2 - r_obs_2^2;
%     Sk_next = (st_next(1)-x_obs_2)^2 + (st_next(2)-y_obs_2)^2 + (st_next(3)-z_obs_2)^2 - (r_obs_2+1)^2 ;
%     temp2 = hk_next/(hk_next-Sk_next);
%     f_bar2 = if_else(temp2 > 0, exp(-1/temp2), 0)/(if_else(temp2 > 0, exp(-1/temp2), 0) + if_else(1-temp2 > 0, exp(-1/(1-temp2)), 0));
%     Phi_next = f_bar2;
%     V_next = [st_next(1) st_next(2) st_next(3) st_next(4)] *P_lyap*[st_next(1); st_next(2); st_next(3); st_next(4)];
%     rho_next = Phi_next/V_next;
% 
%     f_value = f1(st);
%     f_value_next = f1(st_next);
%     C1 = (f_value_next(1)-f_value(1))/(st_next(1)-st(1)) + (f_value_next(2)-f_value(2))/(st_next(2)-st(2))...
%         + (f_value_next(3)-f_value(3))/(st_next(3)-st(3)) + (f_value_next(4)-f_value(4))/(st_next(4)-st(4))...
%         + (f_value_next(5)-f_value(5))/(st_next(5)-st(5)) + (f_value_next(6)-f_value(6))/(st_next(6)-st(6))...
%         + (f_value_next(7)-f_value(7))/(st_next(7)-st(7)) + (f_value_next(8)-f_value(8))/(st_next(8)-st(8));
% 
%     C1 = C1*rho;
% 
%     g_value = g1(st);
%     g_value_next = g1(st_next);
%     C2 = 0;
%     for i=1:4
%         temp = (g_value_next(1,i) - g_value(1,i)/(st(1)-st_next(1))) + (g_value_next(2,i) - g_value(2,i)/(st(2)-st_next(2)))...
%               + (g_value_next(3,i) - g_value(3,i)/(st(3)-st_next(3))) + (g_value_next(4,i) - g_value(4,i)/(st(4)-st_next(4)))...
%               + (g_value_next(5,i) - g_value(5,i)/(st(5)-st_next(5))) + (g_value_next(6,i) - g_value(6,i)/(st(6)-st_next(6)))...
%               + (g_value_next(7,i) - g_value(7,i)/(st(7)-st_next(7))) + (g_value_next(8,i) - g_value(8,i)/(st(8)-st_next(8)));
%         temp = temp * rho * con(i);
%         C2 = C2 + temp;
%     end
% 
%     g = [g; rho_next - rho + C2 + C1];
% 
% end





% make the decision variable one column  vector
OPT_variables = [reshape(X,n_states*(N+1),1);reshape(U,n_controls*N,1)];
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;
args.lbg(1:n_states*(N+1)) = 0; % equality constraints
args.ubg(1:n_states*(N+1)) = 0; % equality constraints

args.lbg(n_states*(N+1)+1 : n_states*(N+1)+ (0*N)) = 0; % inequality constraints
args.ubg(n_states*(N+1)+1 : n_states*(N+1)+ (0*N)) = inf; % inequality constraints

args.lbx(1:n_states:n_states*(N+1),1) = xmin(1); %state x lower bound
args.ubx(1:n_states:n_states*(N+1),1) = xmax(1); %state x upper bound
args.lbx(2:n_states:n_states*(N+1),1) = xmin(2); %state y lower bound
args.ubx(2:n_states:n_states*(N+1),1) = xmax(2); %state y upper bound

args.lbx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = umin(1); %u lower bound
args.ubx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = umax(1); %u upper bound



%% Simulation
t0 = 0;


xlog(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,n_controls);
X0 = repmat(x0,1,N+1)';

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];

tic
while(norm((x0-xf),2) > 1e-2 && mpciter < time_total / dt)
    args.p   = [x0;xf]; % set the values of the parameters vector
    % initial value of the optimization variables
    args.x0  = [reshape(X0',n_states*(N+1),1);reshape(u0',n_controls*N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    u = reshape(full(sol.x(n_states*(N+1)+1:end))',n_controls,N)'; % get controls only from the solution
    xx1(:,1:n_states,mpciter+1)= reshape(full(sol.x(1:n_states*(N+1)))',n_states,N+1)'; % get solution TRAJECTORY
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    % Apply the control and shift the solution
    [t0, x0, u0] = shift(dt, t0, x0, u,f);
    xlog(:,mpciter+1) = x0;
    X0 = reshape(full(sol.x(1:n_states*(N+1)))',n_states,N+1)'; % get solution TRAJECTORY
    % Shift trajectory to initialize the next step
    X0 = [X0(2:end,:);X0(end,:)];
    mpciter
    mpciter = mpciter + 1;
end
toc


Line_width = 2;
Line_color = 'black';

figure
hold on
subplot(2,1,1);
plot(linspace(0, time_total, length(xlog) ), xlog(1,:),...
    'LineWidth',Line_width ,'MarkerSize',4,'Color',Line_color);
% xlabel('$t(s)$','interpreter','latex','FontSize',20);
ylabel('$x(m) $','interpreter','latex','FontSize',10);

subplot(2,1,2);
hold on
plot(linspace(0, time_total, length(xlog(1,:)) ), xlog(2,:),...
    'LineWidth', Line_width,'MarkerSize',4,'Color',Line_color);
% xlabel('$t(s)$','interpreter','latex','FontSize',20);
ylabel('$y(m) $','interpreter','latex','FontSize',10);
% 
% subplot(2,2,3);
% hold on
% plot(linspace(0, time_total, length(xlog(1,:)) ), xlog(3,:),...
%     'LineWidth', Line_width,'MarkerSize',4,'Color',Line_color);
% xlabel('$Time(sec.)$','interpreter','latex','FontSize',10);
% ylabel('$z(m) $','interpreter','latex','FontSize',10);
% 
% subplot(2,2,4);
% hold on
% plot(linspace(0, time_total, length(xlog(1,:)) ), xlog(4,:),...
%     'LineWidth', Line_width,'MarkerSize',4,'Color',Line_color);
% xlabel('$Time(sec.)$','interpreter','latex','FontSize',10);
% ylabel('$\psi(deg)$','interpreter','latex','FontSize',10);
% 
% figure
% hold on
% subplot(2,2,1);
% hold on
% plot(linspace(0, time_total, length(xlog(1,:)) ), xlog(5,:),...
%     'LineWidth', Line_width,'MarkerSize',4,'Color',Line_color);
% % xlabel('$t(s)$','interpreter','latex','FontSize',20);
% ylabel('$u(m/s)$','interpreter','latex','FontSize',10);
% 
% subplot(2,2,2);
% hold on
% plot(linspace(0, time_total, length(xlog(1,:)) ), xlog(6,:),...
%     'LineWidth', Line_width,'MarkerSize',4,'Color',Line_color);
% % xlabel('$t(s)$','interpreter','latex','FontSize',20);
% ylabel('$v(m/s)$','interpreter','latex','FontSize',10);
% 
% subplot(2,2,3);
% hold on
% plot(linspace(0, time_total, length(xlog(1,:)) ), xlog(7,:),...
%     'LineWidth', Line_width,'MarkerSize',4,'Color',Line_color);
% xlabel('$Time(sec.)$','interpreter','latex','FontSize',10);
% ylabel('$w(m/s)$','interpreter','latex','FontSize',10);
% 
% subplot(2,2,4);
% hold on
% plot(linspace(0, time_total, length(xlog(1,:)) ), xlog(8,:),...
%     'LineWidth', Line_width,'MarkerSize',4,'Color',Line_color);
% xlabel('$Time(sec.)$','interpreter','latex','FontSize',10);
% ylabel('$r(rad/s)$','interpreter','latex','FontSize',10);
% % 
% % 
% %    sgtitle('States')
% 
% 
% figure
% hold on
% subplot(2,2,1);
% plot(linspace(0, time_total, length(u_cl(:,1))),u_cl(:,1),...
%     'LineWidth', Line_width,'MarkerSize',4,'Color',Line_color);
% xlabel('$t (s)$','interpreter','latex','FontSize',20);
% ylabel('$f_{surge}(N)$','interpreter','latex','FontSize',10);
% 
% subplot(2,2,2);
% hold on
% plot(linspace(0, time_total, length(u_cl(:,1))),u_cl(:,2),...
%     'LineWidth', Line_width,'MarkerSize',4,'Color',Line_color);
% % xlabel('$t (s)$','interpreter','latex','FontSize',20);
% ylabel('$f_{sway}(N)$','interpreter','latex','FontSize',10);
% 
% subplot(2,2,3);
% hold on
% plot(linspace(0, time_total, length(u_cl(:,1))),u_cl(:,3),...
%     'LineWidth', Line_width,'MarkerSize',4,'Color',Line_color);
% xlabel('$Time(sec.)$','interpreter','latex','FontSize',10);
% ylabel('$f_{heave}(N)$','interpreter','latex','FontSize',10);
% 
% subplot(2,2,4);
% hold on
% plot(linspace(0, time_total, length(u_cl(:,1))),u_cl(:,4),...
%     'LineWidth', Line_width,'MarkerSize',4,'Color',Line_color);
% xlabel('$Time(sec.)$','interpreter','latex','FontSize',10);
% ylabel('$\tau_{yaw}(N.m)$','interpreter','latex','FontSize',10);
% % sgtitle('Control Inputs')
% 
% 
% figure
% % hold on
% plot3(xlog(1,:), xlog(2,:), xlog(3,:),'LineWidth', Line_width,'Color','red')
% xlabel('x(m)','interpreter','latex','FontSize',20);
% ylabel('y(m)','interpreter','latex','FontSize',20);
% zlabel('z(m)','interpreter','latex','FontSize',20);
% hold on
% [X,Y,Z] = sphere;
% % 
% surf(r_obs_1*(X)+x_obs_1,r_obs_1*(Y)+y_obs_1,r_obs_1*(Z)+z_obs_1)
% surf(r_obs_2*(X)+x_obs_2,r_obs_2*(Y)+y_obs_2,r_obs_2*(Z)+z_obs_2)