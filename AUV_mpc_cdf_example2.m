clc;
clear;
close all;
import casadi.*
addpath dynamics\ density_functions\ barrier_functions\

% setup colors for plots
colors = colororder;
blue = colors(1,:);
red = colors(2,:);
yellow = colors(3,:);
green = colors(5,:);
obsColor = [.7 .7 .7]; % Obstacle color -> Grey

%% Setup Parameters
% ---------- system setup ---------------
% states for AUV
x = SX.sym('x');
y = SX.sym('y');
z = SX.sym('z');
psi = SX.sym('psi');
xdot = SX.sym('xdot');
ydot = SX.sym('ydot');
zdot = SX.sym('zdot');
psidot = SX.sym('psidot');
states = [x; y; z; psi; xdot; ydot; zdot; psidot];
n_states = length(states);

% control Inputs for AUV
u = SX.sym('u');
v = SX.sym('v');
w = SX.sym('w');
r = SX.sym('r');
controls = [u;v;w;r];
n_controls = length(controls);

%---------- MPC setup ----------------------
time_total = 10; % time for the total steps, equal to tsim
dt = 0.02;
o = 3;
Q = 10*diag([10,10,10, 10, o, o, o, o]);
P_weight = 100*diag([10,10,10,10, o, o, o, o]);
R = 1*diag([1, 1, 1, 1]);
N = 10;
C_t = 0.1;

xmin = [-inf; -inf; -inf;-100; -100; -100; -100; -100];
xmax = -xmin;

umin = [-inf; -inf; -inf; -inf];
umax = -umin;

% ----------- Environment setup --------------------
x0 = [0; -10; -8; 0; 0; 0; 0; 0]; % initial
x_ini = x0;
xf = [4; 4; 8; 0; 0; 0; 0; 0]; % target

obs_x = SX.sym('obs_x');
obs_y = SX.sym('obs_y');
obs_z = SX.sym('obs_z');
obs_r = SX.sym('obs_r');
obs_s = SX.sym('obs_s');
obs = [obs_x;obs_y;obs_z;obs_r;obs_s];

% obstacle list for complex 3D world
num_obs = 4; % Number of obstacles
obs_toros_xy = [0;0;2;2;3];
obs_toros_xz = [8;0;0;1;2];
obs_cylinder = [-3;0;0;2;3];
obs_sphere = [0;-6;-4;3;10];

% ------------ Density function setup ------------
rho_toros_xy = density_sphere(states,obs_toros_xy);
rho_toros_xy = Function('rho',{states,obs},{rho_toros_xy});

rho_toros_xz = density_sphere(states,obs_toros_xz);
rho_toros_xz = Function('rho',{states,obs},{rho_toros_xz}); 

rho_cylinder = density_sphere(states,obs_cylinder);
rho_cylinder = Function('rho',{states,obs},{rho_cylinder}); 

rho_sphere = density_sphere(states,obs_sphere);
rho_sphere = Function('rho',{states,obs},{rho_sphere}); 



%% Dynamics Setup 
[dx_dt,f,g] = AUV_dynamics(states, controls, dt);
f_discrete = dt*f + states;
g_discrete = dt*g;

% define matlab functions for divergence of f_discrete
jacob_f_discrete = jacobian(f_discrete, states');
div_f_discrete = sum(diag(jacob_f_discrete));
div_f_discrete = Function('Div_F',{states},{div_f_discrete});

% define matlab functions for divergence of g_discrete for each column
g_discrete1 = g_discrete(:,1);
jacob_G_discrete1 = jacobian(g_discrete1, states');
div_g_discrete1 = sum(diag(jacob_G_discrete1));

g_discrete2 = g_discrete(:,2);
jacob_g_discrete2 = jacobian(g_discrete2, states');
div_g_discrete2 = sum(diag(jacob_g_discrete2));

g_discrete3 = g_discrete(:,3);
jacob_g_discrete3 = jacobian(g_discrete3, states');
div_g_discrete3 = sum(diag(jacob_g_discrete3));

g_discrete4 = g_discrete(:,4);
jacob_g_discrete4 = jacobian(g_discrete4, states');
div_g_discrete4 = sum(diag(jacob_g_discrete4));

% calculate the total divergence of g_discrete
div_g_discrete = div_g_discrete1+div_g_discrete2+div_g_discrete3+div_g_discrete4;
div_g_discrete = Function('div_G_discrete1',{states},{div_g_discrete});


% define matlab functions for F=f+gu, f, g
F = Function('F',{states,controls},{dx_dt}); 
f = Function('f',{states},{f}); 
g = Function('g',{states},{g});

%% Casdai MPC setup
% A vector that represents the states over the optimization problem.
X = SX.sym('X',n_states,(N+1));

% Decision variables for control
U = SX.sym('U',n_controls,N); 

% Decision variables for slack
C = SX.sym('C',N); 

% parameters (which include at the initial state of the robot and the reference state)
P = SX.sym('P',n_states + n_states);

obj = 0; % Objective function
constraints = [];  % constraints vector

st  = X(:,1); % initial state
constraints = [constraints;st-P(1:8)]; % initial condition constraints

%------------- Compute cost and constriants -------------------------
% compute running cost and dynamics constraint
for k = 1:N
    st = X(:,k);
    con = U(:,k);
    obj = obj+(st-P(9:16))'*Q*(st-P(9:16)) + con'*R*con; % calculate obj
    obj = obj + sqrt((C(k))^2);
    st_next = X(:,k+1);
    f_value = F(st,con);
    st_next_euler = st+ (dt*f_value);
    constraints = [constraints;st_next-st_next_euler]; % compute constraints
end

% Add Terminal Cost
k = N+1;
st = X(:,k);
obj = obj+(st-P(9:16))'*P_weight*(st-P(9:16)); % calculate obj

% density constraint for obstacles
% add toros in xy
for k = 1:N
    % get current and next state
    st = X(:,k); st_next = X(:,k+1);
    % get current control
    con = U(:,k);

    % get obs location
    obs_loc = obs_toros_xy;

    % get current and next rho
    rho = rho_toros_xy(st,obs_loc);
    rho_next = rho_toros_xy(st_next,obs_loc);
    
    % form density constraint
    density_constraint = (rho_next-rho) + dt*(div_f_discrete(st)+div_g_discrete(st))*rho;
    slack = dt*C(k)*rho;
    constraints = [constraints; density_constraint - slack];
end

% add toros in xz
for k = 1:N
    % get current and next state
    st = X(:,k); st_next = X(:,k+1);
    % get current control
    con = U(:,k);

    % get obs location
    obs_loc = obs_toros_xz;

    % get current and next rho
    rho = rho_toros_xz(st,obs_loc);
    rho_next = rho_toros_xz(st_next,obs_loc);
    
    % form density constraint
    density_constraint = (rho_next-rho) + dt*(div_f_discrete(st)+div_g_discrete(st))*rho;
    slack = dt*C(k)*rho;
    constraints = [constraints; density_constraint - slack];
end

% add cylinder
for k = 1:N
    % get current and next state
    st = X(:,k); st_next = X(:,k+1);
    % get current control
    con = U(:,k);

    % get obs location
    obs_loc = obs_cylinder;

    % get current and next rho
    rho = rho_cylinder(st,obs_loc);
    rho_next = rho_cylinder(st_next,obs_loc);
    
    % form density constraint
    density_constraint = (rho_next-rho) + dt*(div_f_discrete(st)+div_g_discrete(st))*rho;
    slack = dt*C(k)*rho;
    constraints = [constraints; density_constraint - slack];
end

% add sphere
for k = 1:N
    % get current and next state
    st = X(:,k); st_next = X(:,k+1);
    % get current control
    con = U(:,k);

    % get obs location
    obs_loc = obs_sphere;

    % get current and next rho
    rho = rho_sphere(st,obs_loc);
    rho_next = rho_sphere(st_next,obs_loc);
    
    % form density constraint
    density_constraint = (rho_next-rho) + dt*(div_f_discrete(st)+div_g_discrete(st))*rho;
    slack = dt*C(k)*rho;
    constraints = [constraints; density_constraint - slack];
end


%------------- Setup optimization problem -------------------------
% make the decision variable one column  vector
OPT_variables = [reshape(X,n_states*(N+1),1);reshape(U,n_controls*N,1);reshape(C,N,1)];
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', constraints, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;
args.lbg(1:n_states*(N+1)) = 0; % equality constraints
args.ubg(1:n_states*(N+1)) = 0; % equality constraints

args.lbg(n_states*(N+1)+1 : n_states*(N+1)+ (num_obs*N)) = 0; % inequality constraints
args.ubg(n_states*(N+1)+1 : n_states*(N+1)+ (num_obs*N)) = inf; % inequality constraints

args.lbx(1:n_states:n_states*(N+1),1) = xmin(1); %state x lower bound
args.ubx(1:n_states:n_states*(N+1),1) = xmax(1); %state x upper bound
args.lbx(2:n_states:n_states*(N+1),1) = xmin(2); %state y lower bound
args.ubx(2:n_states:n_states*(N+1),1) = xmax(2); %state y upper bound
args.lbx(3:n_states:n_states*(N+1),1) = xmin(3); %state z lower bound
args.ubx(3:n_states:n_states*(N+1),1) = xmax(3); %state z upper bound

args.lbx(4:n_states:n_states*(N+1),1) = xmin(4); %state psi lower bound
args.ubx(4:n_states:n_states*(N+1),1) = xmax(4); %state psi upper bound
args.lbx(5:n_states:n_states*(N+1),1) = xmin(5); %state xdot lower bound
args.ubx(5:n_states:n_states*(N+1),1) = xmax(5); %state xdot upper bound
args.lbx(6:n_states:n_states*(N+1),1) = xmin(6); %state ydot lower bound
args.ubx(6:n_states:n_states*(N+1),1) = xmax(6); %state ydot upper bound
args.lbx(7:n_states:n_states*(N+1),1) = xmin(7); %state zdot lower bound
args.ubx(7:n_states:n_states*(N+1),1) = xmax(7); %state zdot upper bound
args.lbx(8:n_states:n_states*(N+1),1) = xmin(8); %state psidot lower bound
args.ubx(8:n_states:n_states*(N+1),1) = xmax(8); %state psidot upper bound

args.lbx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = umin(1); %u lower bound
args.ubx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = umax(1); %u upper bound
args.lbx(n_states*(N+1)+2:n_controls:n_states*(N+1)+n_controls*N,1) = umin(2); %p lower bound
args.ubx(n_states*(N+1)+2:n_controls:n_states*(N+1)+n_controls*N,1) = umax(2); %p upper bound
args.lbx(n_states*(N+1)+3:n_controls:n_states*(N+1)+n_controls*N,1) = umin(3); %q lower bound
args.ubx(n_states*(N+1)+3:n_controls:n_states*(N+1)+n_controls*N,1) = umax(3); %q upper bound
args.lbx(n_states*(N+1)+4:n_controls:n_states*(N+1)+n_controls*N,1) = umin(4); %r lower bound
args.ubx(n_states*(N+1)+4:n_controls:n_states*(N+1)+n_controls*N,1) = umax(4); %r upper bound

args.lbx(n_states*(N+1)+n_controls*N+1:1:n_states*(N+1)+n_controls*N+N,1) = 0; %C lower bound
args.ubx(n_states*(N+1)+n_controls*N+1:1:n_states*(N+1)+n_controls*N+N,1) = inf; %C upper bound


%% Simulate MPC controller with AUV dynamics
t0 = 0;
xlog(:,1) = x0; % xx contains the history of states
t(1) = t0;
u0 = zeros(N,n_controls);
X0 = repmat(x0,1,N+1)';
C_0 = repmat(C_t,1,N);

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];
C_log = [];

w_bar = waitbar(0,'1','Name','Simulating MPC-CDF...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

while(norm((x0-xf),2) > 1e-2 && mpciter < time_total / dt)
    max_iter = time_total/dt;
    waitbar(mpciter/max_iter,w_bar,sprintf(string(mpciter)+'/'+string(max_iter)))
    
    args.p   = [x0;xf]; % set the values of the parameters vector

    % initial value of the optimization variables
    args.x0  = [reshape(X0',n_states*(N+1),1);reshape(u0',n_controls*N,1);reshape(C_0',N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    u = reshape(full(sol.x(n_states*(N+1)+1:end-10))',n_controls,N)'; % get controls only from the solution
    xx1(:,1:n_states,mpciter+1)= reshape(full(sol.x(1:n_states*(N+1)))',n_states,N+1)'; % get solution TRAJECTORY
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;

    % Apply the control and shift the solution
    [t0, x0, u0] = shift(dt, t0, x0, u,F);
    xlog(:,mpciter+1) = x0;
    X0 = reshape(full(sol.x(1:n_states*(N+1)))',n_states,N+1)'; % get solution TRAJECTORY
    C_0 = reshape(full(sol.x(end-10+1:end))',1,N);
    C_log = [C_log; C_0];

    % Shift trajectory to initialize the next step
    X0 = [X0(2:end,:);X0(end,:)];
    mpciter = mpciter + 1;
end

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% ---------------- plot 3D trajectory ----------------------  

figure(1)

% plot start, goal
plot3(x_ini(1), x_ini(2), x_ini(3), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'black', 'MarkerEdgeColor','black'); hold on;
plot3(xf(1), xf(2), xf(3), 'o', 'MarkerSize', 10, 'MarkerFaceColor', green, 'MarkerEdgeColor', green); hold on;

% plot obstacles
% plot torus in xy
cx = obs_toros_xy(1); cy = obs_toros_xy(2); cz = obs_toros_xy(3); % center
[theta,phi] = meshgrid(linspace(0,2*pi,25));
R = 10; r = obs_toros_xy(4);
x = (R + r*cos(theta)).*cos(phi) + cx;
y = (R + r*cos(theta)).*sin(phi) + cy;
z = r*sin(theta) + cz;
surf(x,y,z); hold on
axis equal

% plot torus in xz
R = 7; r = obs_toros_xz(4);
cx = obs_toros_xz(1); cy = obs_toros_xz(2); cz = obs_toros_xz(3);
x = (R + r*cos(theta)).*cos(phi) + cx;
z = (R + r*cos(theta)).*sin(phi) + cz;
y = r*sin(theta) + cy;
surf(x,y,z); hold on
axis equal

%  plot cylinder
r = obs_cylinder(4); h = 40;
[X,Y,Z] = cylinder(r,20);
Z = Z*h;
surf(X+obs_cylinder(1),Y+obs_cylinder(2),Z-20); hold on
axis equal

% plot shere
[X,Y,Z] = sphere;
r = obs_sphere(4);
X = X*r;
Y = Y*r;
Z = Z*r;
surf(X+obs_sphere(1),Y+obs_sphere(2),Z+obs_sphere(3)); hold on
axis equal

colormap gray

% plot x-y-z trajecotry
plot3(xlog(1,:), xlog(2,:), xlog(3,:),'LineWidth', 2,'Color',red)
hold on

% Set the remaining axes properties
view([-23 16]);
% view([45,25]);
axes1 = gca;
box(axes1,'on');
axis(axes1,'square');
hold(axes1,'off');
grid on;
set(axes1,'FontSize',15,'LineWidth',1.5)
xlabel('$x_1$','interpreter','latex', 'FontSize', 20);
ylabel('$x_2$','interpreter','latex', 'FontSize', 20);
zlabel('$x_3$','interpreter','latex', 'FontSize', 20);

%% ------ plot density function over time
figure
rho_log = full(rho_sphere(xlog,obs_sphere));
time_log = linspace(0, time_total, length(xlog(1,:)));
plot(time_log, rho_log,'LineWidth', 2);
xlabel('$t(s)$','interpreter','latex','FontSize',20);
ylabel('$\rho(x) $','interpreter','latex','FontSize',10);


%% -------------- plots vs time ---------------------------
% % plot states vs time
% figure(2)
% hold on
% subplot(2,2,1);
% plot(linspace(0, time_total, length(xlog) ), xlog(1,:),...
%     'LineWidth',Line_width ,'MarkerSize',4,'Color',Line_color);
% % xlabel('$t(s)$','interpreter','latex','FontSize',20);
% ylabel('$x(m) $','interpreter','latex','FontSize',10);
% 
% subplot(2,2,2);
% hold on
% plot(linspace(0, time_total, length(xlog(1,:)) ), xlog(2,:),...
%     'LineWidth', Line_width,'MarkerSize',4,'Color',Line_color);
% % xlabel('$t(s)$','interpreter','latex','FontSize',20);
% ylabel('$y(m) $','interpreter','latex','FontSize',10);
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
% % plots controls1 vs time
% figure(3)
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
% 
% % controls2 vs time
% figure(4)
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
