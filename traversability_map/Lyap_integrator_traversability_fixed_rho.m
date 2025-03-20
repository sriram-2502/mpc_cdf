clc;
clear;
close all;
import casadi.*
addpath ..\dynamics\ ..\density_functions\ ..\barrier_functions\ ..\utils\

% Set defaults
set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultAxesColor', 'w');
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesLineWidth', 1);

% setup colors for plots
colors = colororder;
blue = colors(1,:);
red = colors(2,:);
yellow = colors(3,:);
green = colors(5,:);
obsColor = [.7 .7 .7]; % Obstacle color -> Grey

%% Setup Parameters
% ---------- system setup -------------------------------------------------
% states for integrator
x = SX.sym('x');
y = SX.sym('y');

states = [x; y];
n_states = length(states);

% control Inputs in density space
u1 = SX.sym('u1');
u2 = SX.sym('u2');
controls = [u1;u2];
n_controls = length(controls);

%---------- MPC setup ----------------------
time_total = 10; % time for the total steps, equal to tsim
tracking = 0; % set to 1 to track a ref traj
N = 10;
dt = 0.1; dt_sim = 0.1;
Q = 10*diag([1,1]);
R = 1*diag([1, 1]);
P_terminal = 10*diag([1,1]); % terminal cost
P_trav = 1e2; % weight on trav cost
C_t = 0.1;

xmin = [-inf; -inf];
xmax = -xmin;
umin = [-1; -1];
umax = -umin;

% ------------- env setup -------------------------------------------------
% initial Conditions on a grid
x0 = [1;13]; x_ini = x0;
xf = [25;13]; % target

% ------------- load height map and fit rbfs ------------------------------
grid_map = readmatrix('Terrain Map.xlsx','Sheet','hill_and_pit');
max_val = max(grid_map,[],"all");
height_map = grid_map./max_val;
[rows, cols] = size(height_map);

% Generate the uniform grid centers
grid_spacing = 1; sigma = 1; 
[center_x, center_y] = meshgrid(1:grid_spacing:cols, 1:grid_spacing:rows);
centers = [center_x(:), center_y(:)];
% Fit RBFs and get the weights
weights = fit_rbf(height_map, centers, sigma);

% ------------ Density function setup ------------
b_x = query_rbf_height(states', weights, centers, sigma)+5;
V_x = b_x;
lyap_fun = Function('rho',{states},{V_x});

%% Dynamics Setup (for divergence constraint)
% dynamics without paramter mismatch
[dx_dt,f,g] = integrator_dynamics(states, controls);
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

% calculate the total divergence of g_discrete
div_g_discrete = div_g_discrete1+div_g_discrete2;%+div_g_discrete3+div_g_discrete4;
div_g_discrete = Function('div_G_discrete1',{states},{div_g_discrete});

% define matlab functions for F=f+gu, f, g
F = Function('F',{states,controls},{dx_dt}); 
f = Function('f',{states},{f}); 
g = Function('g',{states},{g});

%% Casdai MPC setup
% A vector that represents the states over the optimization problem.
X = SX.sym('X',n_states,(N+1));

% Decision variables for control (rho bar)
U = SX.sym('U',n_controls, N); 

% Decision variables for slack
C = SX.sym('C',N); 

% parameters (which include at the initial state of the robot and the reference state)
if(tracking)
    % for ref tracking
    P = SX.sym('P',n_states + N*n_states); 
else
    % for point stabliziation
    P = SX.sym('P',n_states + n_states); 
end

obj = 0; % Objective function
constraints = [];  % constraints vector

% CONSTRAINT: initial condition constraints
st  = X(:,1); % initial state
constraints = [constraints;st-P(1:n_states)]; 

%------------- Compute cost and constriants -------------------------
% compute running cost and dynamics constraint
for k = 1:N
    st = X(:,k);
    con = U(:,k);
    
    % COST: get cost function
    % sum{q(x)*rho + rho_bar.T*R*rho_bar/rho}
    if(tracking)
        % for tracking
        q_x = (st-P(n_states*k+1:n_states*k+n_states))'*Q*(st-P(n_states*k+1:n_states*k+n_states));
    else
        % for point stabilization
        q_x = (st-P(n_states+1:2*n_states))'*Q*(st-P(n_states+1:2*n_states));
    end
%     obj = obj + q_x*rho + (con'*R*con)/rho;
    obj = obj + q_x + (con'*R*con);
    
    % CONSTRAINT: get dynamics constraint
    % x(k+1) = F(x(k),u(k))
    st_next = X(:,k+1);
    f_value = F(st,con);
    st_next_euler = st + (dt*f_value);
    constraints = [constraints;st_next-st_next_euler]; 
end

% Add Terminal Cost
if(~tracking)  
    k = N+1;
    st = X(:,k);
    terminal_cost = (st-P(n_states+1:2*n_states))'*P_terminal*(st-P(n_states+1:2*n_states)); % calculate obj
    obj = obj + terminal_cost;
end

% CONSTRAINT: Delta(V) < 0 for stability
for k = 1:N
    % get current and next state
    st = X(:,k); st_next = X(:,k+1);
    
    % get current and next rho
    V = lyap_fun(st');
    V_next = lyap_fun(st_next');

    % form density constraint
    alpha = 1;
    lyap_stability = V_next - V + alpha*V;
    constraints = [constraints; lyap_stability];
end

% CONSTRAINT/Obj:  b(x_k)*rho_k <= gamma
traversability = 0;
num_centers = size(centers, 1);
rbf_values = SX.zeros(1, num_centers);
for k = 1:N
    st = X(:,k); 
    for i = 1:num_centers
        diff = st' - centers(i, :);
        rbf_values(i) = exp(-sum(diff.^2, 2) / (2 * sigma^2));
    end
    b_x = rbf_values * weights; % height estimate
%     traversability = traversability + b_x*rho;
    traversability = traversability + b_x;
end
% constraints = [constraints; traversability];
obj = obj + P_trav*traversability;

%------------- Setup optimization problem -------------------------
% make the decision variable one column  vector
OPT_variables = [reshape(X,n_states*(N+1),1); reshape(U,n_controls*N,1);...
    reshape(C,N,1);];
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', constraints, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);
args = struct;

%------------------- equality constraints for dyanmics -------------------------------
args.lbg(1:n_states*(N+1)) = 0; 
args.ubg(1:n_states*(N+1)) = 0;

%------------------ inequality constraints -------------------------------
% bounds for CONSTRAINT: divergence constriant from MPC-CDF paper
args.lbg(n_states*(N+1)+1 : n_states*(N+1)+N) = -Inf; 
args.ubg(n_states*(N+1)+1 : n_states*(N+1)+N) = -1e-2; 

% bounds for CONSTRAINT: b(x(k))*rho(k) <= gamma
% one constraint for entire horizon (sum over horizon)
% args.lbg(n_states*(N+1)+N+1) = 0; 
% args.ubg(n_states*(N+1)+N+1) = gamma; 

% ------------------- bounds ----------------------------------------------
args.lbx(1:n_states:n_states*(N+1),1) = xmin(1); %state x lower bound
args.ubx(1:n_states:n_states*(N+1),1) = xmax(1); %state x upper bound
args.lbx(2:n_states:n_states*(N+1),1) = xmin(2); %state y lower bound
args.ubx(2:n_states:n_states*(N+1),1) = xmax(2); %state y upper bound

args.lbx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = umin(1); %u1 lower bound
args.ubx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = umax(1); %u1 upper bound
args.lbx(n_states*(N+1)+2:n_controls:n_states*(N+1)+n_controls*N,1) = umin(2); %u2 lower bound
args.ubx(n_states*(N+1)+2:n_controls:n_states*(N+1)+n_controls*N,1) = umax(2); %u2 upper bound

args.lbx(n_states*(N+1)+n_controls*N+1:n_states*(N+1)+n_controls*N+N,1) = 0; %C lower bound
args.ubx(n_states*(N+1)+n_controls*N+1:n_states*(N+1)+n_controls*N+N,1) = inf; %C upper bound

%% Simulate MPC controller
% init
t0 = 0;
t(1) = t0;
u0 = zeros(N,n_controls);
X0 = repmat(x0,1,N+1)';
C_0 = repmat(C_t,1,N);
tstart = 0; tend = dt_sim;

% logging
x_sol = [];
u_cl=[];
rho_sol = [];
C_log = [];
tlog = [];
xlog(:,1) = x0; 
V_log = [];
V_log = [V_log, lyap_fun(x0')];
tlog(1) = t0;

% Start MPC
mpciter = 1;
w_bar = waitbar(0,'1','Name','Simulating MPC-CDF...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

while(norm((x0-xf),2) > 1e-2 && mpciter < time_total / dt_sim)
    max_iter = time_total/dt_sim;
    current_time = mpciter*dt_sim;
    waitbar(mpciter/max_iter,w_bar,sprintf(string(mpciter)+'/'+string(max_iter)))
    
    % set the values of the parameters vector
    if(tracking)
        % get ref trjaectory for states
        args.p(1:n_states) = x0;
%         for k=1:N
%             t_predict = current_time + (k-1)*dt_sim;
%             x_ref = x0(1) + 2*t_predict;
%             y_ref = x0(2);
%             if(x_ref >= 25)
%                 xref  = 25;
%             end
%             args.p(n_states*k+1:n_states*k+n_states) = [x_ref, y_ref];
%         end
        t_predict = current_time + N*dt_sim;
        x_target = [2*t_predict; 13];
        x_ref = generate_reference(x0, x_target, N, dt);
        args.p(n_states+1:n_states+N*n_states) = x_ref(:);
    else
        args.p   = [x0;xf]; % for point stabilization
    end
    
    % initial value of the optimization variables
    args.x0  = [reshape(X0',n_states*(N+1),1);reshape(u0',n_controls*N,1);...
        reshape(C_0',N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p);
    
    % get solutions
    u = reshape(full(sol.x(n_states*(N+1)+1:n_states*(N+1)+n_controls*N))',n_controls,N)'; 
    C_0 = reshape(full(sol.x(end-N+1:end))',1,N);   
    
    u_star = u(1,:);
    [t,X] = ode45(@(t,x)full(F(x,u_star)),[tstart,tend],x0);
    % --- update ---
    x0 = X(end,:)';
    tstart = tend; 
    tend = tend + dt_sim;
    t0 = tend;
    
    % logging
    xlog(:,mpciter+1) = x0;
    V_log = [V_log, lyap_fun(x0)];
    u_cl= [u_cl ; u(1,:)];
    tlog(mpciter+1) = t0;
    C_log = [C_log; C_0];
    
    % Shift trajectory to initialize the next step
    X0 = reshape(full(sol.x(1:n_states*(N+1)))',n_states,N+1)'; 
    X0 = [X0(2:end,:);X0(end,:)];
    mpciter = mpciter + 1;
end

waitbar_graphics = findall(0,'type','figure','tag','TMWWaitbar');
delete(waitbar_graphics);

%% ---------------- plot 2D trajectory ----------------------  
close all
figure(1)
subplot(1,2,1)
% plot heightmap
[x, y] = meshgrid(1:cols, 1:rows);
surf(x, y, height_map, 'FaceAlpha', 0.5);
xlabel('X');
ylabel('Y');
zlabel('Height');
view(2)
hold on

% plot x-y-z trajecotry
traj = plot(xlog(1,:), xlog(2,:),'-','LineWidth', 2,'Color',red);
xlabel('x(m)','interpreter','latex','FontSize',20);
ylabel('y(m)','interpreter','latex','FontSize',20);
hold on

% plot start and target 
plot(x_ini(1), x_ini(2), 'o', 'MarkerSize',10, 'MarkerFaceColor','black','MarkerEdgeColor','black'); hold on;
if(~tracking)
    plot(xf(1), xf(2), 'o', 'MarkerSize',10, 'MarkerFaceColor',green,'MarkerEdgeColor',green); hold on;
end
%setup plots
axes1 = gca;
box(axes1,'on');
axis(axes1,'square');

hold(axes1,'off');
xlabel('position, $x$ (m)','interpreter','latex', 'FontSize', 20);
ylabel('position, $y$ (m)','interpreter','latex', 'FontSize', 20);
xlim([1,28]); ylim([1,28])

%%%%%%%%%%%%%%
subplot(1,2,2)
% plot heightmap reconstructed using rbfs
% Generate a fine grid of coordinates
[x_fine, y_fine] = meshgrid(linspace(1, cols, 10*cols), linspace(1, rows, 10*rows));
coordinates_fine = [x_fine(:), y_fine(:)];

% Reconstruct the height map using the RBFs and weights on the fine grid
num_centers = size(centers, 1);
rbf_matrix_fine = zeros(size(coordinates_fine, 1), num_centers);
for i = 1:num_centers
    diff = coordinates_fine - centers(i, :);
    rbf_matrix_fine(:, i) = exp(-sum(diff.^2, 2) / (2 * sigma^2));
end
fitted_height_map_fine = reshape(rbf_matrix_fine * weights, size(x_fine));

surf(x_fine, y_fine, fitted_height_map_fine, 'FaceAlpha', 0.5, 'EdgeColor','none');
xlabel('X');
ylabel('Y');
zlabel('Height');
view(2)
hold on

% plot x-y-z trajecotry
traj = plot(xlog(1,:), xlog(2,:),'-','LineWidth', 2,'Color',red);
xlabel('x(m)','interpreter','latex','FontSize',20);
ylabel('y(m)','interpreter','latex','FontSize',20);
hold on

% plot start and target 
plot(x_ini(1), x_ini(2), 'o', 'MarkerSize',10, 'MarkerFaceColor','black','MarkerEdgeColor','black'); hold on;
if(~tracking)
    plot(xf(1), xf(2), 'o', 'MarkerSize',10, 'MarkerFaceColor',green,'MarkerEdgeColor',green); hold on;
end

%setup plots
axes1 = gca;
box(axes1,'on');
axis(axes1,'square');
xlim([1,28]); ylim([1,28])
hold(axes1,'off');
xlabel('position, $x$ (m)','interpreter','latex', 'FontSize', 20);
ylabel('position, $y$ (m)','interpreter','latex', 'FontSize', 20);

%% plot ctrol and density functions vs time
figure(2)
subplot(3,2,[1,2])
plot(tlog, xlog(1,:),'-','LineWidth', 2, 'Color',blue); hold on;
plot(tlog, xlog(2,:),'LineWidth', 2,'Color',red);
xlabel('time (s)','interpreter','latex', 'FontSize', 20);
ylabel('states (m)','interpreter','latex', 'FontSize', 20);
lgd = legend('x','y');
lgd.Interpreter = 'latex';
lgd.FontSize = 15;
grid on
xlim([0,time_total])

subplot(3,2,[3,4])
plot(tlog(1:end-1), u_cl(:,1),'LineWidth', 2,'Color',blue); hold on;
plot(tlog(1:end-1), u_cl(:,2),'LineWidth', 2,'Color',red);
xlabel('time (s)','interpreter','latex', 'FontSize', 20);
ylabel('control (m/s)','interpreter','latex', 'FontSize', 20);
lgd = legend('$u_1$','$u_2$');
lgd.Interpreter = 'latex';
lgd.FontSize = 15;
grid on
xlim([0,time_total])

subplot(3,2,[5,6])
plot(tlog, full(V_log), 'LineWidth', 2,'Color',blue);
xlabel('time (s)','interpreter','latex', 'FontSize', 20);
ylabel('Lyapunov, $V(x)$','interpreter','latex', 'FontSize', 20);
lgd = legend('$V$');
lgd.Interpreter = 'latex';
lgd.FontSize = 15;
grid on
xlim([0,time_total])