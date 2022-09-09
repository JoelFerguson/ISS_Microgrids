function dxdt = DC_pH(t,x)
global n m B R L C qd Gl Il Kp A Kd1 E Ki Kd3 P1



P = diag(P1);


% Physical states
phi = x(1:n,1);
q = x(n+1:2*n,1);
%phit = x(2*n+1:2*n+m,1);

X = zeros(n,n);
Gleq = zeros(n,n);
Delta = zeros(n,1);


for i=1:n
        X(i,i) = inv(qd(i)/C(i,i))*P(i,i)*inv(q(i)/C(i,i));
        Gleq(i,i) = Gl(i,i) - X(i,i);
        Delta(i) = Il(i) + inv(qd(i)/C(i,i))*P(i,i) + Gl(i,i)*qd(i)/C(i,i);
end


% Controller state
xc = x(2*n+1:3*n,1);


% Auxiliary control input   
v = inv(C)*(q-qd) + R*inv(L)*phi - Kd1*Kp*(phi+A*(q-qd))...
     -E*A*Kp*(phi+A*(q-qd)) - E*inv(C)*(q-qd)...
     -(Kd1+E*A-0.5*A*inv(L)*inv(Kp))*Ki*(phi+A*(q-qd)-xc);


% Control input
u = inv(C)*qd + v;


% Closed loop dynamics

dxdt(1:n,1) = - R*inv(L)*phi - inv(C)*q + u; %phi

dxdt(n+1:2*n,1) = inv(L)*phi - Gleq*inv(C)*(q-qd)...
                     - Delta; %q
                
%dxdt(2*n+1:2*n+m,1) = - (B')*inv(C)*q - Rt*inv(Lt)*phit; % phit

dxdt(2*n+1:3*n,1) = -0.5*inv(Kp)*inv(L)*A*Kp*(phi+A*(q-qd))...
                        -inv(Kp)*inv(L)*inv(C)*(q-qd)...
                        +(Kd3-0.5*inv(Kp)*inv(L)*A)*Ki*(phi+A*(q-qd)-xc); % xc




