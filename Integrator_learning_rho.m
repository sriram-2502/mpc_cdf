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
% states for integrator
x = SX.sym('x');
y = SX.sym('y');
rho = SX.sym('rho');
states = [x; y;];
n_states = length(states);

% control Inputs in density space
rho_bar1 = SX.sym('rho_bar1');
rho_bar2 = SX.sym('rho_bar2');
controls = [rho_bar1;rho_bar2];
n_controls = length(controls);

%---------- MPC setup ----------------------
time_total = 20; % time for the total steps, equal to tsim
N = 10; % for mismatch use N = 100
dt = 0.01; % use dt = 0.1 for cbf and vanilla obs
o = 1;
Q = 10*diag([10,10]);
R = 1*diag([1, 1]);
C_t = 0.1;

xmin = [-inf; -inf];
xmax = -xmin;

umin = [-inf; -inf];
umax = -umin;

% % ----------- Environment setup --------------------
% % initial Conditions on a grid
% x0 = [0;0.01]; x_ini = x0;
% xf = [10;0]; % target
% 
% obs_x = SX.sym('obs_x');
% obs_y = SX.sym('obs_y');
% obs_r = SX.sym('obs_r');
% obs_s = SX.sym('obs_s');
% obs = [obs_x;obs_y;obs_r;obs_s];
% 
% % obstacle list for sphere world
% num_obs = 1; % Number of obstacles
% obs_rad = 1;
% obs_sens = obs_rad + 1;
% obs1 = [4; 0; obs_rad(1); obs_sens(1)];
% 
% %---------- cbf/obstalce constraint setup ---------------
% rho_circle = density_circle(states,obs);
% rho_circle = Function('b',{states,obs},{rho_circle}); 



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

% Decision variables for density function values
rho = SX.sym('rho',1,(N+1));

% Decision variables for control (rho bar)
rho_bar = SX.sym('U',n_controls, N); 

% Decision variables for slack
C = SX.sym('C',N); 

% parameters (which include at the initial state of the robot and the reference state)
% P = SX.sym('P',n_states + n_states); % for point stabliziation
P = SX.sym('P',n_states + N*n_states); % for ref tracking

obj = 0; % Objective function
constraints = [];  % constraints vector

st  = X(:,1); % initial state
constraints = [constraints;st-P(1:n_states)]; % initial condition constraints

%------------- Compute cost and constriants -------------------------
% compute running cost and dynamics constraint
for k = 1:N
    st = X(:,k);
    con = U(:,k);
%     obj = obj+(st-P(n_states+1:2*n_states))'*Q*(st-P(n_states+1:2*n_states)) + con'*R*con; % calculate obj
    obj = obj+(st-P(n_states*k+1:n_states*k+n_states))'*Q*(st-P(n_states*k+1:n_states*k+n_states)) + con'*R*con;
    
    st_next = X(:,k+1);
    f_value = F(st,con);
    st_next_euler = st+ (dt*f_value);
    constraints = [constraints;st_next-st_next_euler]; % x(k+1) = F(x(k),u(k))
end

% Add Terminal Cost
k = N+1;
st = X(:,k);
% density constraint for obstacles
for obs_num = 1:num_obs
    for k = 1:N
        % get current and next state
        st = X(:,k); st_next = X(:,k+1);
        % get current control
        con = U(:,k);

        % get obs location
        obs_loc = eval(sprintf('obs%d',obs_num));

        % get current and next rho
        rho = rho_circle(st,obs_loc);
        rho_next = rho_circle(st_next,obs_loc);
        
        % form density constraint
        density_constraint = (rho_next-rho) + dt*(div_f_discrete(st)+div_g_discrete(st))*rho;
        slack = dt*C(k)*rho;
        constraints = [constraints; density_constraint - slack];
    end
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

args.lbx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = umin(1); %u1 lower bound
args.ubx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = umax(1); %u1 upper bound
args.lbx(n_states*(N+1)+2:n_controls:n_states*(N+1)+n_controls*N,1) = umin(2); %u2 lower bound
args.ubx(n_states*(N+1)+2:n_controls:n_states*(N+1)+n_controls*N,1) = umax(2); %u2 upper bound

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
    current_time = mpciter*dt;
    waitbar(mpciter/max_iter,w_bar,sprintf(string(mpciter)+'/'+string(max_iter)))
    
    % set the values of the parameters vector
%     args.p   = [x0;xf]; 
    args.p(1:n_states) = x0;
    for k=1:N
        t_predict = current_time + (k-1)*dt;
        x_ref = 0.5*t_predict;
        y_ref = 1;
        if(x_ref >=12)
            xref  = 12;
            y_ref = 12;
        end
        args.p(n_states*k+1:n_states*k+n_states) = [x_ref, y_ref];
    end    
    
    % initial value of the optimization variables
    args.x0  = [reshape(X0',n_states*(N+1),1);reshape(u0',n_controls*N,1);reshape(C_0',N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p);
    u = reshape(full(sol.x(n_states*(N+1)+1:end-N))',n_controls,N)'; % get controls only from the solution
    xx1(:,1:n_states,mpciter+1)= reshape(full(sol.x(1:n_states*(N+1)))',n_states,N+1)'; % get solution TRAJECTORY
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    
    C_0 = reshape(full(sol.x(end-N+1:end))',1,N);
    C_log = [C_log; C_0];

    % Apply the control and shift the solution
    [t0, x0, u0] = shift(dt, t0, x0, u,F);

    xlog(:,mpciter+1) = x0;
    X0 = reshape(full(sol.x(1:n_states*(N+1)))',n_states,N+1)'; % get solution TRAJECTORY

    % Shift trajectory to initialize the next step
    X0 = [X0(2:end,:);X0(end,:)];
    mpciter = mpciter + 1;
end

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% ---------------- plot 2D trajectory ----------------------  

figure(1)
% For legend as rectangular object can't be defined as a legend
    dummy_marker = plot(NaN,NaN, 'o','MarkerSize', 10, 'MarkerEdgeColor',...
            'black', 'MarkerFaceColor',obsColor, 'LineWidth', 1.5); 
        
% plot obstacles
xc = obs1(1); yc = obs1(2); Rc = obs1(3);
angles = (0:100-1)*(2*pi/100);
points = [xc;yc] + [Rc*cos(angles);Rc*sin(angles)];
P = polyshape(points(1,:), points(2,:));
plot(P, 'FaceColor', obsColor, 'LineWidth', 2, 'FaceAlpha', 1.0); hold on;

% plot x-y-z trajecotry
traj = plot(xlog(1,:), xlog(2,:),'LineWidth', 2,'Color',red);
xlabel('x(m)','interpreter','latex','FontSize',20);
ylabel('y(m)','interpreter','latex','FontSize',20);
hold on

% plot start and target 
% plot(x_ini(1), x_ini(2), 'o', 'MarkerSize',10, 'MarkerFaceColor','black','MarkerEdgeColor','black'); hold on;
% plot(xf(1), xf(2), 'o', 'MarkerSize',10, 'MarkerFaceColor',green,'MarkerEdgeColor',green); hold on;


%setup plots
axes1 = gca;
box(axes1,'on');
axis(axes1,'equal');

% set other axis properties
% set(axes1,'FontSize',15,'LineWidth',2);
% if(cbf_constraint)
%     lgd = legend(traj,'CBF with $\gamma$ = '+string(gamma));
% else
%     lgd = legend(traj,'distance');
% end

lgd.Interpreter = 'latex';
lgd.FontSize = 15;
grid on;
xlim([0,10])
ylim([yc-5,yc+5])
hold(axes1,'off');
xlabel('Position, $x$ (m)','interpreter','latex', 'FontSize', 20);
ylabel('Position, $y$ (m)','interpreter','latex', 'FontSize', 20);