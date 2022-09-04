%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: DCnetwork_idealSim.m
% Description: Microgrid simulation example for "Increasing the region of
% attraction in DC microgrids" using ideal network dynamics
% Authour: Michele Cucuzzella, Joel Ferguson
% Date 2-September-2022
% Version: 1.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

% Function to remove NaNs when starting from 0 voltage
removeNan = @(x) fillmissing(x,'constant',0);

%% Simulation settings
% Select simulation scenario
% 1: Ideal power loads
% 2: Non-ideal power loads
simNumber = 2; 

% Simulation length
sim.Tf = 9;
% Time of power load step change
sim.Tload1 = 3;
% Time to add time-varying power load
sim.Tload2 = 6;

%% Define plant dynamcis
% n nodes, m edges
n = 4;
m = 4;

% Network incidence matrix
sys.B = [-1 0 0 -1;       
      1 -1 0 0;
      0  1 -1 0;
    0 0 1 1];

% Define filter and network parameters
sys.R = diag([0.2,0.3,0.5,0.1]);
sys.L = diag([1.8,2,3,2.2])*1e-3;
sys.C = diag([2.2,1.9,2.5,1.7])*1e-3;
sys.Rt = diag([70,50,80,60])*1e-3;
sys.Lt = diag([2.1,2.3,2,1.8])*1e-6;

% Define nominal load parameters
sys.GLs = diag([1/10,1/6,1/4,1/8]);
sys.ILs = [4 4 4 4]';
sys.PLs = @(t) [250 250 250 250]' + (t > sim.Tload1)*[750 0 0 0]';
% Define load cut-off voltages for two simulated scenarios
switch simNumber
    case 1 % Ideal constant power load; Vc = 0
        sys.Vc = zeros(n,1);% cut-off voltage (ideal CPL: sys.Vc = 0, non ideal CPL: sys.Vc = 24)
        sys.GLp = @(t) zeros(n);
        sim.V0 = [23.95; 24; 23.9; 23.85];

    case 2 % Non-ideal constant power load; Vc = 24
        sys.Vc = 24*ones(n,1);
        sys.GLp = @(t) diag(sys.PLs(t)./sys.Vc.^2);
        sim.V0 = zeros(n,1);
end
% Define ZIP load based on defined parameters
sys.IL = @(q,t) removeNan((sys.C\q >= sys.Vc).*(sys.GLs*(sys.C\q) + sys.ILs + sys.PLs(t)./(sys.C\q))) + (sys.C\q < sys.Vc).*((sys.GLs + sys.GLp(t))*(sys.C\q) + sys.ILs);

% Define time-varing disturbance which is unknown to the controller
sys.delta = @(q,t) zeros(2*n+m,1) + removeNan((t > sim.Tload2)*[zeros(n,1); -100*sin(t)*sys.C(1,1)/q(1); zeros(n-1,1); zeros(m,1)]);

% Define functinos to partition state vector for x = [phi; q; phi_t] 
sys.phi = @(x) x(1:n);
sys.q = @(x) x(n+1:2*n);
sys.phit = @(x) x(2*n+1:2*n+m);

% Define open-loop plant dynamics
sys.H = @(phi,q,phit) 0.5*phi.'*(sys.L\phi) + 0.5*q.'*(sys.C\q) + 0.5*phit.'*(sys.Lt\phit);
sys.dHdx = @(phi,q,phit) [sys.L\phi; sys.C\q; sys.Lt\phit];
sys.F = [-sys.R -eye(n) zeros(n,m);
            eye(n) zeros(n) sys.B;
            zeros(m,n) -sys.B.' -sys.Rt];
sys.dx = @(t,x,u) sys.F*sys.dHdx(sys.phi(x),sys.q(x),sys.phit(x)) + [u; -sys.IL(sys.q(x),t); zeros(m,1)] + sys.delta(sys.q(x),t);

%% Define controller
% Target note voltages and charges
ctrl.V_star = [47.9 48 47.8 47.7].';
ctrl.q_star = sys.C*ctrl.V_star;

% Define feedforward controller from proposition 2
ctrl.u = @(v) sys.C\ctrl.q_star + v;

% Define parameters for ISS controller
% Tuning parameters
ctrl.Kp = 4e1*eye(n);
ctrl.Ki = 1e0*eye(n);
ctrl.A = 1.5e0*eye(n);
ctrl.Kd1 = 1.5e2*ctrl.A;
% Computed parameters
ctrl.E = inv(sys.L)*inv(ctrl.Kp) + inv(sys.L)*ctrl.A*sys.C*ctrl.A;
ctrl.Kd2 = inv(sys.L)*ctrl.A*sys.C;
ctrl.Kd3 = 0.25*inv(ctrl.Kp)*inv(sys.L)*ctrl.A*inv(ctrl.Kd1+ctrl.E*ctrl.A)*ctrl.A*inv(sys.L)*inv(ctrl.Kp);

% Gains for implementation
% Integrator gains
ctrl.w_x = (ctrl.Kd3 - 0.5*inv(ctrl.Kp)*inv(sys.L)*ctrl.A)*ctrl.Ki;
ctrl.w_q = 0.5*inv(sys.L)*ctrl.A*ctrl.A + inv(ctrl.Kp)*inv(sys.L)*inv(sys.C) - ctrl.w_x*ctrl.A;
ctrl.w_phi = 0.5*inv(sys.L)*ctrl.A - ctrl.w_x;
% Signal gains
ctrl.psi_x = -(ctrl.Kd1 + ctrl.E*ctrl.A - 0.5*ctrl.A*inv(sys.L)*inv(ctrl.Kp))*ctrl.Ki;
ctrl.psi_q = -inv(sys.C) + ctrl.Kd1*ctrl.Kp*ctrl.A + ctrl.E*ctrl.Kp*ctrl.A*ctrl.A + ctrl.E*inv(sys.C) - ctrl.psi_x*ctrl.A;
ctrl.psi_phi = -sys.R*inv(sys.L) + ctrl.Kd1*ctrl.Kp + ctrl.E*ctrl.A*ctrl.Kp - ctrl.psi_x;

% Define ISS controller from Proposition 3
ctrl.dxc = @(q_tilde,phi,xc) -ctrl.w_q*q_tilde - ctrl.w_phi*phi - ctrl.w_x*xc;
ctrl.v = @(q_tilde,phi,xc) -ctrl.psi_q*q_tilde - ctrl.psi_phi*phi - ctrl.psi_x*xc;

% Partition function to extract controller state
ctrl.xc = @(x) x(2*n+m+1:3*n+m);
% Function to compute the charge error
ctrl.q_tilde = @(x) x(n+1:2*n) - ctrl.q_star;

% Rewrite controller to expect full state vector as input
ctrl.dxc_wrap = @(x) ctrl.dxc(ctrl.q_tilde(x), sys.phi(x), ctrl.xc(x));
ctrl.v_wrap = @(x) ctrl.v(ctrl.q_tilde(x), sys.phi(x), ctrl.xc(x));

%% Construct ISS-Lyapunov function for closed-loop system
% Define lumped terms from Proposition 3
ctrl.Phi = -(sys.Rt\sys.B.')*(sys.C\ctrl.q_star);
ctrl.Lambda = @(q,t) (sys.C\q >= sys.Vc).*(sys.ILs + sys.PLs(t)./(sys.C\ctrl.q_star) + sys.GLs*(sys.C\ctrl.q_star) - sys.B*ctrl.Phi) + (sys.C\q < sys.Vc).*(sys.ILs + (sys.GLs + sys.GLp(t))*(sys.C\ctrl.q_star) - sys.B*ctrl.Phi);
ctrl.Gamma = @(q,t) ctrl.Ki\inv(ctrl.Kd1 + ctrl.E*ctrl.A - 0.5*ctrl.A*inv(sys.L)*inv(ctrl.Kp))*(ctrl.Kd1 + ctrl.E*ctrl.A)*ctrl.Kp*sys.L*ctrl.Lambda(q,t);

% Compute the equilibrium values of lumped terms
ctrl.Lambda_star = @(t) ctrl.Lambda(ctrl.q_star,t);
ctrl.Gamma_star = @(t) ctrl.Gamma(ctrl.q_star,t);

% Compute the equilibrium point for given load and input
ctrl.X_star = @(t) [sys.L*ctrl.Lambda_star(t); ctrl.q_star; sys.Lt*ctrl.Phi; sys.L*ctrl.Lambda_star(t)+ctrl.Gamma_star(t)];


% Define the ISS-Lyapunov function
ctrl.Q = [ctrl.Kp+ctrl.Ki, (ctrl.Kp+ctrl.Ki)*ctrl.A, zeros(n,m), -ctrl.Ki;
        ctrl.A*(ctrl.Kp+ctrl.Ki), ctrl.A*(ctrl.Kp+ctrl.Ki)*ctrl.A+inv(sys.C) zeros(n,m), -ctrl.A*ctrl.Ki;
        zeros(m,n), zeros(m,n), inv(sys.Lt), zeros(m,n);
        -ctrl.Ki, ctrl.Ki*ctrl.A, zeros(n,m), ctrl.Ki];
ctrl.W = @(X,t) 0.5*(X - ctrl.X_star(t)).'*ctrl.Q*(X - ctrl.X_star(t));

%% Simulate system
% Define initial conditions
sim.phi0 = zeros(n,1);
sim.phit0 = zeros(m,1);
sim.xc0 = zeros(n,1);
sim.x0 = [sim.phi0; (sys.C*sim.V0); sim.phit0; sim.xc0];

% Pack ODE for closed-loop system
sim.wrap = @(t,x) [sys.dx(t,x,ctrl.u(ctrl.v_wrap(x)));
                    ctrl.dxc_wrap(x)];

% Simulate system
options = odeset('RelTol',1e-6);
[res.t, res.x] = ode23s(sim.wrap,[0 sim.Tf],sim.x0,options);

% Compute result
res.V = zeros(length(res.t),n);
for i=1:length(res.t)
    q = res.x(i,n+1:2*n).';
    res.W(i) = ctrl.W(res.x(i,:).',res.t(i));
    res.V(i,:) = (sys.C\q).';
end

%% Plot results
% Node voltages
figure(1)
plot(res.t, res.V, 'LineWidth', 2)
grid on
legend('V_1','V_2','V_3','V_4')

% Lyapunov function
figure(2)
plot(res.t, res.W, 'LineWidth', 2)
grid on
legend('Lyapunov function')
