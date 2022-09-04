%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: DCnetwork_nonIdealSim.m
% Description: Microgrid simulation example for "Increasing the region of
% attraction in DC microgrids" using non-ideal network dynamics using the Simscape Electrical toolbox to model switching behaviours. 
% Authour: Michele Cucuzzella, Joel Ferguson
% Date 4-September-2022
% Version: 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

%% simulation settings
flag1 = 1; % 0: constant load, 1: time varying load
Vc = 1*24*ones(4,1);% cut-off voltage (ideal CPL: Vc = 0, non ideal CPL: Vc = 24)

Tf = 9; % simulation time
Tload1 = 3; % time instant when P changes
Tload2 = 6; % time instant when P changes
Dtl = 0.01; % P takes Dtl seconds for changing from P1 to P2 

% Parameters for converting voltage reference into PWM signals
Ts = 1e-5;%1e-5
f = 1e4;%1e4
upper_limit_sat = 0.99;

%% Network parameters
% n nodes, m edges
n = 4;
m = 4;

% Parameters for switches in voltage source
Ron = 1e-6;
Rs = 1e6;
% DC voltage source level
Vdc = diag([100 100 100 100]);

% incidence matrix
B = [-1 0 0 -1;       
      1 -1 0 0;
      0  1 -1 0;
    0 0 1 1];

% filter and network parameters
 R = diag([0.2,0.3,0.5,0.1]);
 L = diag([1.8,2,3,2.2])*1e-3;
 C = diag([2.2,1.9,2.5,1.7])*1e-3;
Rt = diag([70,50,80,60])*1e-3;
Lt = diag([2.1,2.3,2,1.8])*1e-6;

% load parameters
Gl = diag([1/10,1/6,1/4,1/8]);
Glp = zeros(n,n);
Il = [4 4 4 4]';
P1 = [250 250 250 250];
P2 = (P1 + [3 0 0 0]*250);
% Time varying load parameters
sine_freq = 5;
if flag1 == 1
    sine_ampl = 0.1e3;
else
    sine_ampl = 0;
end

% desired voltages
Vd = [47.9 48 47.8 47.7];
% desired charges
qd = C*Vd';

% controller parameters
Kp = 5e1*eye(n);%5e1
Ki = 1e0*eye(n);%1e0
A = 1.5e0*eye(n);%1.5e0
Kd1 = 1.5e2*A;%1.5e2
E = inv(L)*inv(Kp) + inv(L)*A*C*A;
Kd3 = 0.25*inv(Kp)*inv(L)*A*inv(Kd1+E*A)*A*inv(L)*inv(Kp);
Kd2 = inv(L)*A*C;

%% Check that parametwers and initial conditions satisfy requirements
% % check on the equivalent conductance at the steady-state
% check1 = eig(Gl - inv(diag(Vd)^2)*diag(P1))
% check2 = eig(Gl - inv(diag(Vd)^2)*diag(P2))

% % check on the controller parameters (considering P2 is the worst case)
% Gleq = zeros(n,n);
% for i=1:n
%     if Vc(i)== 0
%         Gleq(i,i) = Gl(i,i) - inv(qd(i)/C(i,i))*P2(i)*inv(q0(i)/C(i,i));
%     else
%         Gleq(i,i) = Gl(i,i) - inv(qd(i)/C(i,i))*P2(i)/Vc(i);
%     end
% end
% check3 = eig(A + 2*L*Gleq*inv(C)) % positive definite (32a)

% if 1.5*A*Gl*inv(2*Gleq+Kd2)*Gl*A >= 1.5*A*B*inv(Rt)*B'*A
%     check4 = eig(Kd1-1.5*A*Gl*inv(2*Gleq+Kd2)*Gl*A) % positive definite (32b)
% else
%     check4 = eig(Kd1-1.5*A*B*inv(Rt)*B'*A) % positive definite (32b)
% end

%% Define simulation intiial conditions
% initial voltage conditions
if Vc > 0
    V0 = 0*Vd;
else
    V0 = 0.5*Vd;
end
q0 = C*V0';

% flux initial conditions
Glp0 = zeros(n,n);
phitbar = - Lt*inv(Rt)*B'*inv(C)*q0;
for i=1:n
    if Vc(i)>0 % possible resistive behavior
        Glp0(i,i) = P1(1)/(Vc(i))^2;
    end
    if q0(i)>=C(i,i)*Vc(i) % CPL behavior
        phibar = L*(- B*inv(Lt)*phitbar + Gl*inv(C)*q0 + Il + inv(diag(inv(C)*q0))*P1');
    else % resistive behavior
        phibar = L*(- B*inv(Lt)*phitbar + (Gl+Glp0)*inv(C)*q0 + Il);
    end
end

% Integrator initial conditions
xcbar = -0.5*inv(Ki)*inv(Kd3-0.5*inv(Kp)*inv(L)*A)*inv(Kp)*inv(L)*A*Kp*phibar + phibar;

%% Simulate system
sim('DCnetworkModel')

%% Plot results
% Intial transient response 0-0.1s
figure(1)
subplot(2,2,1)
[~,endIdx] = min(abs(V.time-0.1));
plot(V.time(1:endIdx),V.signals.values(1:endIdx,1),'LineWidth',2)
hold on
plot(V.time(1:endIdx),Vd(1)*ones(size(V.time(1:endIdx))),'--','LineWidth',2)
grid on
ylabel('V_1 [V]')

subplot(2,2,2)
[~,endIdx] = min(abs(V.time-0.1));
plot(V.time(1:endIdx),V.signals.values(1:endIdx,2),'LineWidth',2)
hold on
plot(V.time(1:endIdx),Vd(2)*ones(size(V.time(1:endIdx))),'--','LineWidth',2)
grid on
ylabel('V_2 [V]')

subplot(2,2,3)
[~,endIdx] = min(abs(V.time-0.1));
plot(V.time(1:endIdx),V.signals.values(1:endIdx,3),'LineWidth',2)
hold on
plot(V.time(1:endIdx),Vd(3)*ones(size(V.time(1:endIdx))),'--','LineWidth',2)
grid on
ylabel('V_3 [V]')

subplot(2,2,4)
[~,endIdx] = min(abs(V.time-0.1));
plot(V.time(1:endIdx),V.signals.values(1:endIdx,4),'LineWidth',2)
hold on
plot(V.time(1:endIdx),Vd(4)*ones(size(V.time(1:endIdx))),'--','LineWidth',2)
grid on
ylabel('V_4 [V]')


% Response to step-change and time-varyng loads
figure(2)
subplot(2,2,1)
[~,startIdx] = min(abs(V.time-2));
plot(V.time(startIdx:end),V.signals.values(startIdx:end,1),'LineWidth',2)
hold on
plot(V.time(startIdx:end),Vd(1)*ones(size(V.time(startIdx:end))),'--','LineWidth',2)
xlim([2 9])
grid on
ylabel('V_1 [V]')

subplot(2,2,2)
[~,startIdx] = min(abs(V.time-2));
plot(V.time(startIdx:end),V.signals.values(startIdx:end,2),'LineWidth',2)
hold on
plot(V.time(startIdx:end),Vd(2)*ones(size(V.time(startIdx:end))),'--','LineWidth',2)
xlim([2 9])
grid on
ylabel('V_2 [V]')

subplot(2,2,3)
[~,startIdx] = min(abs(V.time-2));
plot(V.time(startIdx:end),V.signals.values(startIdx:end,3),'LineWidth',2)
hold on
plot(V.time(startIdx:end),Vd(3)*ones(size(V.time(startIdx:end))),'--','LineWidth',2)
xlim([2 9])
grid on
ylabel('V_3 [V]')

subplot(2,2,4)
[~,startIdx] = min(abs(V.time-2));
plot(V.time(startIdx:end),V.signals.values(startIdx:end,4),'LineWidth',2)
hold on
plot(V.time(startIdx:end),Vd(4)*ones(size(V.time(startIdx:end))),'--','LineWidth',2)
xlim([2 9])
grid on
ylabel('V_4 [V]')




