% Automatica v1

clc
clear all
close all


% global variables
global n m B R L C Vd qd Gl Il Kp A Kd1 E Ki Kd3 P1

% n nodes, m edges
n = 1;
m = 0;

Tf = 5; % simulation time

% incidence matrix
B = 0;

% desired voltages
Vd = 48;

% filter and network parameters
 R = 0.2;
 L = 1.8*1e-3;
 C = 2.2*1e-3;
Lt = 0;

% desired charges
qd = C*Vd';

q0 = qd;

% load parameters
Gl = 0;
Il = 0;
P1 = 5*250;

% check on the equivalent conductance at the steady-state
check1 = eig(Gl - inv(diag(Vd)^2)*diag(P1))

% controller parameters
Kp = 4e1*eye(n);
Ki = 1e0*eye(n);
A = 0.9e0*eye(n);%0.9;0.9e1
Kd1 = 0.01e0*A;%0.01
E = inv(L)*inv(Kp) + inv(L)*A*C*A;
Kd3 = 0.25*inv(Kp)*inv(L)*A*inv(Kd1+E*A)*A*inv(L)*inv(Kp);
Kd2 = inv(L)*A*C;

% check on the controller parameters
Gleq = zeros(n,n);
for i=1:n
    Gleq(i,i) = Gl(i,i) - inv(qd(i)/C(i,i))*P1(i)*inv(q0(i)/C(i,i));
end

check3 = eig(A + 2*L*Gleq*inv(C)) % positive definite (32a)


check4 = eig(Kd1-1.5*A*Gl*inv(2*Gleq+Kd2)*Gl*A) % positive definite (32b)


%% small RoA phi1, q1 
close all

figure(1)
        ylabel('$q$','FontSize', 18, 'Interpreter', 'latex'); figure(gcf)
        xlabel('$\phi$', 'FontSize', 18, 'Interpreter', 'latex');
        grid on
        hold on
        ay = gca;
ay.TickLabelInterpreter = 'latex';
set(ay,'FontSize', 14)

eps = 0.01;
q0min = q0 - 0.9999*q0;
q0max = 2*q0;
num_sim = 75;
size_mesh = (q0max-q0min)/num_sim;

ylim([0, .3])
xlim([-.4, 1])

%plot(qd(1),qd(2),'ro','Linewidth',2,'MarkerSize',10) % qd

for i = q0min:size_mesh:q0max
        Q0 = i;
        % suitable initial condition        
        phibar = L*(Gl*inv(C)*Q0 + Il + inv(diag(inv(C)*Q0))*P1');
        xcbar = -0.5*inv(Ki)*inv(Kd3-0.5*inv(Kp)*inv(L)*A)*inv(Kp)*inv(L)*A*Kp*phibar + phibar;
        
        x0=[phibar' ... % phi
            Q0'   ... % q
            xcbar'];

        options = odeset('RelTol',1e-3, 'AbsTol',1e-3);    
        [t,x]=ode23t(@(t,x)DC_pH_1node(t,x),[0 Tf],x0,options);

        if abs(x(end,2)-qd(1)) < eps 
            plot(x(:,1),x(:,2), 'Color',[0.8 0.8 0.8], 'LineWidth',1)
            plot(phibar,Q0,'bx','Linewidth',1,'MarkerSize',10) % Q0
            plot(x(end,1),x(end,2), 'ro','Linewidth',2,'MarkerSize',10) % end point
        else
            plot(x(1:22,1),x(1:22,2), 'k', 'LineWidth',1)
            plot(phibar,Q0,'rx','Linewidth',1,'MarkerSize',10) % Q0
        end
end


%% large RoA phi1, q1 
close all
clc

% controller parameters
Kp = 4e1*eye(n);
Ki = 1e0*eye(n);
A = 0.9e1*eye(n);%0.9;0.9e1
Kd1 = 0.01e0*A;%0.01
E = inv(L)*inv(Kp) + inv(L)*A*C*A;
Kd3 = 0.25*inv(Kp)*inv(L)*A*inv(Kd1+E*A)*A*inv(L)*inv(Kp);
Kd2 = inv(L)*A*C;


figure(2)
        ylabel('$q$','FontSize', 18, 'Interpreter', 'latex'); figure(gcf)
        xlabel('$\phi$', 'FontSize', 18, 'Interpreter', 'latex');
        grid on
        hold on
        ay = gca;
ay.TickLabelInterpreter = 'latex';
set(ay,'FontSize', 14)

eps = 0.01;
q0min = q0 - 0.9999*q0;
q0max = 2*q0;
num_sim = 75;
size_mesh = (q0max-q0min)/num_sim;

ylim([0, .3])
xlim([-1, 2.5])

%plot(qd(1),qd(2),'ro','Linewidth',2,'MarkerSize',10) % qd

for i = q0min:size_mesh:q0max
        Q0 = i;
        % suitable initial condition        
        phibar = L*(Gl*inv(C)*Q0 + Il + inv(diag(inv(C)*Q0))*P1');
        xcbar = -0.5*inv(Ki)*inv(Kd3-0.5*inv(Kp)*inv(L)*A)*inv(Kp)*inv(L)*A*Kp*phibar + phibar;
        
        x0=[phibar' ... % phi
            Q0'   ... % q
            xcbar'];

        options = odeset('RelTol',1e-3, 'AbsTol',1e-3);    
        [t,x]=ode23t(@(t,x)DC_pH_1node(t,x),[0 Tf],x0,options);

        if abs(x(end,2)-qd(1)) < eps 
            plot(x(:,1),x(:,2), 'Color',[0.8 0.8 0.8], 'LineWidth',1)
            plot(phibar,Q0,'bx','Linewidth',1,'MarkerSize',10) % Q0
            plot(x(end,1),x(end,2), 'ro','Linewidth',2,'MarkerSize',10) % end point
        else
            plot(x(1:22,1),x(1:22,2), 'k', 'LineWidth',1)
            plot(phibar,Q0,'rx','Linewidth',1,'MarkerSize',10) % Q0
        end
end

 
 
%% small RoA Voltage zoom and Tf
close all

figure(3)
        ylabel('$q$','FontSize', 18, 'Interpreter', 'latex'); figure(gcf)
        xlabel('time', 'FontSize', 18, 'Interpreter', 'latex');
        grid on
        hold on
        ay = gca;
ay.TickLabelInterpreter = 'latex';
set(ay,'FontSize', 14)

eps = 0.01;
q0min = q0 - 0.9999*q0;
q0max = 2*q0;
num_sim = 100;
size_mesh = (q0max-q0min)/num_sim;

xlim([0, Tf])
ylim([0, .213])

%plot(qd(1),qd(2),'ro','Linewidth',2,'MarkerSize',10) % qd

for i = q0min:size_mesh:q0max
        Q0 = i;
        % suitable initial condition
        %phitbar = - Lt*inv(Rt)*B'*inv(C)*Q0;
        
        phibar = L*(Gl*inv(C)*Q0 + Il + inv(diag(inv(C)*Q0))*P1');
        
        xcbar = -0.5*inv(Ki)*inv(Kd3-0.5*inv(Kp)*inv(L)*A)*inv(Kp)*inv(L)*A*Kp*phibar + phibar;
        
        x0=[phibar' ... % phi
            Q0'   ... % q
            xcbar'];

        options = odeset('RelTol',1e-3, 'AbsTol',1e-3);    
        [t,x]=ode23t(@(t,x)DC_pH_1node(t,x),[0 Tf],x0,options);

        if abs(x(end,2)-qd(1)) < eps 
            %plot(Q0(1),Q0(2),'bx','Linewidth',2,'MarkerSize',10) % Q0
            plot(t,x(:,2), 'b', 'LineWidth',1)
            %plot(x(end,3),x(end,4), 'ro','Linewidth',2,'MarkerSize',10) % end point
        else
            %plot(Q0(1),Q0(2),'rx','Linewidth',2,'MarkerSize',10) % Q0
            plot(t,x(:,2), 'm', 'LineWidth',1)
        end
end

plot(t,qd(1)*ones(1,length(t)),'c', 'LineWidth', 2)


%% large RoA Voltage zoom and Tf
close all
clc

% controller parameters
Kp = 4e1*eye(n);
Ki = 1e0*eye(n);
A = 0.9e1*eye(n);%0.9;0.9e1
Kd1 = 0.01e0*A;%0.01
E = inv(L)*inv(Kp) + inv(L)*A*C*A;
Kd3 = 0.25*inv(Kp)*inv(L)*A*inv(Kd1+E*A)*A*inv(L)*inv(Kp);
Kd2 = inv(L)*A*C;

figure(4)
        ylabel('$q$','FontSize', 18, 'Interpreter', 'latex'); figure(gcf)
        xlabel('time', 'FontSize', 18, 'Interpreter', 'latex');
        grid on
        hold on
        ay = gca;
ay.TickLabelInterpreter = 'latex';
set(ay,'FontSize', 14)

eps = 0.01;
q0min = q0 - 0.9999*q0;
q0max = 2*q0;
num_sim = 75;
size_mesh = (q0max-q0min)/num_sim;

xlim([0, 1.5])
ylim([0.01, .3])

%plot(qd(1),qd(2),'ro','Linewidth',2,'MarkerSize',10) % qd

for i = q0min:size_mesh:q0max
        Q0 = i;
        % suitable initial condition
        %phitbar = - Lt*inv(Rt)*B'*inv(C)*Q0;
        
        phibar = L*(Gl*inv(C)*Q0 + Il + inv(diag(inv(C)*Q0))*P1');
        
        xcbar = -0.5*inv(Ki)*inv(Kd3-0.5*inv(Kp)*inv(L)*A)*inv(Kp)*inv(L)*A*Kp*phibar + phibar;
        
        x0=[phibar' ... % phi
            Q0'   ... % q
            xcbar'];

        options = odeset('RelTol',1e-3, 'AbsTol',1e-3);    
        [t,x]=ode23t(@(t,x)DC_pH_1node(t,x),[0 Tf],x0,options);

        if abs(x(end,2)-qd(1)) < eps 
            %plot(Q0(1),Q0(2),'bx','Linewidth',2,'MarkerSize',10) % Q0
            plot(t,x(:,2), 'b', 'LineWidth',1)
            %plot(x(end,3),x(end,4), 'ro','Linewidth',2,'MarkerSize',10) % end point
        else
            %plot(Q0(1),Q0(2),'rx','Linewidth',2,'MarkerSize',10) % Q0
            plot(t,x(:,2), 'm', 'LineWidth',1)
        end
end

plot(t,qd(1)*ones(1,length(t)),'c', 'LineWidth', 2)






