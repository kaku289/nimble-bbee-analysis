tic;
params.Kp = -20;
params.Ki = 0;
params.Kd = 0;
params.a0 = 0;
params.v0 = 0.7;
params.y0 = -1;
params.tau = 0.05;
params.rref = -0.5;
params.y1_minustau = params.y0;

tspan = [0 2.4];
sol = ddensd(@(t,y,ydel,ypdel) ddefun(t,y,ydel,ypdel,params), ...
             params.tau, params.tau, [params.y0 params.v0]', tspan);

tn = linspace(0,2.4);
yn = deval(sol,tn); 

toc;

plot(yn(1,:), yn(2,:));
xlabel('y_1');
ylabel('y_2');


function yp = ddefun(t,y,ydel,ypdel,params)
% Set of neutral-delay-differential equations governing closed-loop landing
% dynamics of bbees

% y is a state vector (=[y; Vgy])
% yp is dy/dt
% ydel and ypdel contain delays in state and state derivatives
% ydel and ypdel have only one column corresponding to tau (sensory delay)
% params - a custom structure containing values of 9 parameters: Kp, Ki,
% Kd, tau, rref, y1_minustau, a0, v0, y0

y1dot = y(2);
y2dot = params.Kp*(params.rref-ydel(2)/ydel(1)) + ...
        params.Ki*(params.rref*t-log(abs(ydel(1)))+log(abs(params.y1_minustau))) - ...
        params.Kd*(ypdel(2)./ydel(1)-(ydel(2)./ydel(1)).^2);

yp = [y1dot;
      y2dot];
end
