function [track, figureHandles] = estimateLandingParameters(track, createPlots)

%   estimateLandingParameters
%   Learns parameters of the landing model based on the derived system of
%   neutral-delayed-differentional equations

%   track - instance of trackForLandingModel
%   creatPlots - true or false
%   outputSignalForEstimation - 'r' or 'ay' (works only for ay at the
%   moment)


%%

% Put measured data in form suitable for system identification

% Create params structure - It contains variables that DO NOT change for a
% track during optimization (for parallel computations)
params.time = track.state(:,1)-track.state(1,1);
params.iodata = [track.state(:,3), track.state(:,6), track.state(:,9)]; % [y Vgy ay]
params.negTime_state = [track.negTime_state(:,1)-track.negTime_state(end,1) track.negTime_state(:,3) ...
                 track.negTime_state(:,6) track.negTime_state(:,9)]; % [negative_time y Vgy ay]
params.y0 = params.iodata(1,1:2)'; % y(0), y is a state vector (= [y; Vgy])
% params.weights = ones(length(params.time),1); % weights used during cost function

% Finding Gaussian weights
r = params.iodata(:,2)./params.iodata(:,1);
% Find time instant that are within +- 0.5 of rref
constr_times = params.time(abs(track.rref-r) <= 0.5);
% Set middle time point as mu and sigma = mu/2
mu = constr_times(ceil(end/2));
sigma = mu/2;
params.weights = exp(-((params.time-mu).^2./2./sigma.^2).^2); % weights used during cost function

params.cost = 'sae'; % can be sae or sse
params.tspan = [params.time(1) params.time(end)];
params.lb = [-20,-20,-1,0.001,track.rref-1]'; % NOTE: Order of parameters is [Kp, Ki, Kd, tau, rref]
params.ub = [20,20,1,0.06,track.rref+1]';
params.ddeoptions = ddeset('Events',@events); % option set used for integrating system of neutral-delayed-differential equations
params.epsilon = 1e-10; % perturbation in parameters to estimate derivatives

% % Creating history solution
% params.histsol.solver = 'ddensd'; 
% params.histsol.x = params.negTime_state(:,1)';
% params.histsol.y = params.negTime_state(:,2:3)';
% params.histsol.yp = params.negTime_state(:,4)';
% params.histsol.IVP = false;
% params.histsol.history = params.negTime_state(1,2:3)';
% params.histsol.xe = []; params.histsol.ye = []; params.histsol.ie = [];

% keyboard;


%% Testing the integration and cost function process (NOT fool proof yet!)
% % params.time = track.state(:,1)-track.state(1,1);
% % params.iodata = [track.state(:,3), track.state(:,6), track.state(:,9)]; % [y Vgy ay]
% % params.negTime_state = [track.negTime_state(:,1)-track.negTime_state(end,1) track.negTime_state(:,3) ...
% %                  track.negTime_state(:,6) track.negTime_state(:,9)]; % [negative_time y Vgy ay]
% % 
% % params.Kp = -20;
% % params.Ki = 0;
% % params.Kd = 0;
% % 
% % params.tau = 0.05;
% % params.rref = -0.5;
% % 
% % params.y0 = [-1 0.7]'; % y(0), y is a state vector (= [y; Vgy])
% % 
% % params.y_minustau = params.y0; % y(-tau), y is a state vector (= [y; Vgy])
% % 
% % params.weights = ones(length(params.time),1); % weights used during cost function
% % params.cost = 'sae'; % can be sae or sse
% % params.tspan = [0 2.1];
% % 
% % % Calling cost function
% % tic;
% % vals = landingDynamics_Objective(params.Kp, params.Ki, params.Kd, params.tau, params.rref, params);
% % toc;
% % keyboard
%% Running the optimization (Smooth solvers)

% creating problem structure
problem = createOptimProblem('fmincon',...
    'objective',@(x) landingDynamics_Objective(x(1), x(2), x(3), x(4), x(5), params),...
    'x0',[0 0 0 0.025 0]',...
    'lb',params.lb,'ub',params.ub, 'options', ...
    optimoptions(@fmincon,'UseParallel',true,'Algorithm','interior-point','Display','off')); % 

% [x,fval] = fmincon(problem)
tic;
gs = GlobalSearch('Display','iter')
gs.StartPointsToRun = 'bounds';
gs.NumTrialPoints = 2000;
% gs.MaxTime = 300;
% rng(14,'twister') % for reproducibility
[x,fval] = run(gs,problem)
toc;
% 
% keyboard;
% 
% % for estimating two parameters - Kp and tau
% problem = createOptimProblem('fmincon',...
%     'objective',@(x) landingDynamics_Objective(x(1), 0, 0, x(2), track.rref, params),...
%     'x0',[-5 0.025]',...
%     'lb',params.lb([1 4]),'ub',params.ub([1 4]), 'options', ...
%     optimoptions(@fmincon,'UseParallel',true,'Algorithm','interior-point','Display','off')); % 
% 
% % [x,fval] = fmincon(problem)
% tic;
% gs = GlobalSearch('Display','iter')
% gs.StartPointsToRun = 'bounds';
% gs.NumTrialPoints = 1000;
% % gs.MaxTime = 300;
% % rng(14,'twister') % for reproducibility
% [x,fval] = run(gs,problem)
% toc;
% Kp = x(1); Ki = 0; Kd = 0; tau = x(2); rref = track.rref;
% 
% %%%%%%%%%%%%%%%%%%%%%%% For user-defined gradient %%%%%%%%%%%%%%%%
% % creating problem structure
% x = [0.374 -7.4 -0.125 0.0263 -0.0312];
% Kp = x(1); Ki = x(2); Kd = x(3); tau = x(4); rref = x(5);
% cf = landingDynamics_Objective(Kp, Ki, Kd, tau, rref, params)
% 
% 
% problem = createOptimProblem('fmincon',...
%     'objective',@(x) landingDynamics_scalarObjective_Jacobian(x(1), x(2), x(3), x(4), x(5), params),...
%     'x0',[0 0 0 0.25 0]',...
%     'lb',params.lb,'ub',params.ub, 'options', ...
%     optimoptions(@fmincon,'UseParallel',true,'Algorithm','interior-point','Display','off' ...
%                  ,'SpecifyObjectiveGradient',true)); %    
% %                  ,'CheckGradients',true,'SpecifyObjectiveGradient',true,'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-3)); % 
% 
% % [tx,fval] = fmincon(problem)
% tic;
% gs = GlobalSearch('Display','iter')
% gs.StartPointsToRun = 'bounds';
% gs.NumTrialPoints = 5000;
% % gs.MaxTime = 300;
% % rng(14,'twister') % for reproducibility
% [x,fval] = run(gs,problem)
% toc;
% 
% keyboard;
%% Running the optimization (Non-smooth solvers)

%%%%%%%%%%%%%%%%%% Pattern search %%%%%%%%%%%%%%%%%%%%%%%
% creating the problem structure
% options = optimoptions('patternsearch','UseCompletePoll',true,'UseParallel',true,'Display','iter');
% 
% % Estimate 5 parameters
% tic; 
% [x,cost,eflag,output] = patternsearch(@(x) landingDynamics_Objective(x(1), x(2), x(3), x(4), x(5), params),...
%     [0 0 0 0.025 0]',[],[],[],[],params.lb,params.ub,[],options); % [0 0 0 0.025 0] [0 0 0 0.025 -1.8] [0.374 -7.4 -0.125 0.0263 -1.8]
% toc;
% Kp = x(1); Ki = x(2); Kd = x(3); tau = x(4); rref = x(5);
% 
% % Estimate 4 parameters
% % rref = -0.833;
% % tic; 
% % [x,cost,eflag,output] = patternsearch(@(x) landingDynamics_Objective(x(1), x(2), x(3), x(4), rref, params),...
% %     [0.374 -7.4 -0.125 0.0263]',[],[],[],[],params.lb,params.ub,[],options); % [0 0 0 0.025 0] [0 0 0 0.025 -1.8] [0.374 -7.4 -0.125 0.0263 -1.8]
% % toc;
% 
% % Estimate Kp, tau and rref
% Ki = 0;
% Kd = 0;
% tic; 
% [x,cost,eflag,output] = patternsearch(@(x) landingDynamics_Objective(x(1), Ki, Kd, x(2), x(3), params),...
%     [-10 0.001 track.rref]',[],[],[],[],[-10 0.0 -5],[10 0.1 0],[],options); % [0 0 0 0.025 0] [0 0 0 0.025 -1.8] [0.374 -7.4 -0.125 0.0263 -1.8]
% toc;
% Kp = x(1); tau = x(2); rref = x(3);
% 
% % Estimate Kp, Ki, tau and rref
% Kd = 0;
% tic; 
% [x,cost,eflag,output] = patternsearch(@(x) landingDynamics_Objective(x(1), x(2), Kd, x(3), x(4), params),...
%     [0 0 0.01 -5]',[],[],[],[],[-30 -30 0.01 -10],[30 30 0.06 0],[],options); % [0 0 0 0.025 0] [0 0 0 0.025 -1.8] [0.374 -7.4 -0.125 0.0263 -1.8]
% toc;
% Kp = x(1); Ki = x(2); tau = x(3); rref = x(4);
% 
% % Estimate Kp, Ki, Kd and tau
% rref = -4;
% tic; 
% [x,cost,eflag,output] = patternsearch(@(x) landingDynamics_Objective(x(1), x(2), x(3), x(4), rref, params),...
%     [0.374 -7.4 -0.125 0.0263]',[],[],[],[],[-10 -10 -1 0.01],[10 10 -1 0.06],[],options); % [0 0 0 0.025 0] [0 0 0 0.025 -1.8] [0.374 -7.4 -0.125 0.0263 -1.8]
% toc;
% Kp = x(1); Ki = x(2); Kd = x(3); tau = x(4);
% 
% % Estimate Kp and tau
% tic; 
% [x,cost,eflag,output] = patternsearch(@(x) landingDynamics_Objective(x(1), 0, 0, x(2), track.rref, params),...
%     [-10 0.08]',[],[],[],[],[-20 0],[20 0.1],[],options); % [0 0 0 0.025 0] [0 0 0 0.025 -1.8] [0.374 -7.4 -0.125 0.0263 -1.8]
% toc;
% Kp = x(1); Ki = 0; Kd = 0; tau = x(2); rref = track.rref;
% 
% % Estimate Kp, Ki and tau
% tic; 
% [x,cost,eflag,output] = patternsearch(@(x) landingDynamics_Objective(x(1), x(2), 0, x(3), track.rref, params),...
%     [-5 -5 0.025]',[],[],[],[],[-10 -10 0.01],[10 10 0.06],[],options); % [0 0 0 0.025 0] [0 0 0 0.025 -1.8] [0.374 -7.4 -0.125 0.0263 -1.8]
% toc;
% Kp = x(1); Ki = x(2); Kd = 0; tau = x(3); rref = track.rref;

% keyboard;

%% Optimization using non-smooth solvers (specifically pattern search with different initial points)
% % Estimate Kp, tau
options = optimoptions('patternsearch','UseCompletePoll',true,'UseParallel',true,'Display','final','MaxIterations',1000);
Ki = 0;
Kd = 0;
rref = track.rref;

[X1, X2] = ndgrid([-5 0 5],[0.015 0.03 0.045]);
set_x0 = [X1(:) X2(:)];

track.param_estimation = outputForParameterEstimation.empty;
tic;
for ct=1:4%size(set_x0,1)
    
    x0 = set_x0(ct,:)';
    
    [xOpt,cost,eflag,output] = patternsearch(@(x) landingDynamics_Objective(x(1), Ki, Kd, x(2), rref, params),...
        x0,[],[],[],[],params.lb([1 4]),params.ub([1 4]),[],options); % [0 0 0 0.025 0] [0 0 0 0.025 -1.8] [0.374 -7.4 -0.125 0.0263 -1.8]
    
    track.param_estimation(ct) = outputForParameterEstimation();
    track.param_estimation(ct).x0 = x0;
    track.param_estimation(ct).xOpt = xOpt;
    track.param_estimation(ct).cost = cost;
    track.param_estimation(ct).eflag = eflag;
    track.param_estimation(ct).output = output;
end
toc;
Kp = track.param_estimation(ct).xOpt(1); tau = track.param_estimation(ct).xOpt(2);

figureHandles = [];

% % Estimate Kp, tau and rref
% options = optimoptions('patternsearch','UseCompletePoll',true,'UseParallel',true,'Display','final','MaxIterations',1000);
% Ki = 0;
% Kd = 0;
% 
% [X1, X2] = ndgrid([-5 0 5],[0.015 0.03 0.045]);
% set_x0 = [X1(:) X2(:)];
% set_x0(:,end+1) = track.rref*ones(size(set_x0,1),1);
% 
% track.param_estimation = outputForParameterEstimation.empty;
% tic;
% for ct=1:4%size(set_x0,1)
%     
%     x0 = set_x0(ct,:)';
%     
%     [xOpt,cost,eflag,output] = patternsearch(@(x) landingDynamics_Objective(x(1), Ki, Kd, x(2), x(3), params),...
%         x0,[],[],[],[],params.lb([1 4 5]),params.ub([1 4 5]),[],options); % [0 0 0 0.025 0] [0 0 0 0.025 -1.8] [0.374 -7.4 -0.125 0.0263 -1.8]
%     
%     track.param_estimation(ct) = outputForParameterEstimation();
%     track.param_estimation(ct).x0 = x0;
%     track.param_estimation(ct).xOpt = xOpt;
%     track.param_estimation(ct).cost = cost;
%     track.param_estimation(ct).eflag = eflag;
%     track.param_estimation(ct).output = output;
% end
% toc;
% Kp = track.param_estimation(ct).xOpt(1); tau = track.param_estimation(ct).xOpt(2); rref = track.param_estimation(ct).xOpt(3);
% 
% figureHandles = [];


% % Estimate Kp, Ki tau and rref
% options = optimoptions('patternsearch','UseCompletePoll',true,'UseParallel',true,'Display','iter','MaxIterations',1000);
% Kd = 0;
% 
% [X1, X2] = ndgrid([-5 0 5],[0.015 0.03 0.045]);
% set_x0 = [X1(:) X2(:)];
% set_x0(:,end+1) = track.rref*ones(size(set_x0,1),1);
% set_x0 = [-2.5*ones(size(set_x0,1),1) set_x0];
% 
% track.param_estimation = outputForParameterEstimation.empty;
% tic;
% for ct=1:4%size(set_x0,1)
%     
%     x0 = set_x0(ct,:)';
%     
%     [xOpt,cost,eflag,output] = patternsearch(@(x) landingDynamics_Objective(x(1), x(2), Kd, x(3), x(4), params),...
%         x0,[],[],[],[],params.lb([1 2 4 5]),params.ub([1 2 4 5]),[],options); % [0 0 0 0.025 0] [0 0 0 0.025 -1.8] [0.374 -7.4 -0.125 0.0263 -1.8]
%     
%     track.param_estimation(ct) = outputForParameterEstimation();
%     track.param_estimation(ct).x0 = x0;
%     track.param_estimation(ct).xOpt = xOpt;
%     track.param_estimation(ct).cost = cost;
%     track.param_estimation(ct).eflag = eflag;
%     track.param_estimation(ct).output = output;
% end
% toc;
% Kp = track.param_estimation(ct).xOpt(1); Ki = track.param_estimation(ct).xOpt(2); tau = track.param_estimation(ct).xOpt(3); rref = track.param_estimation(ct).xOpt(4);
% 
% figureHandles = [];


% % Estimate Kp, Ki, Kd, tau and rref
% options = optimoptions('patternsearch','UseCompletePoll',true,'UseParallel',true,'Display','iter','MaxIterations',1000);
% 
% [X1, X2] = ndgrid([-5 0 5],[0.015 0.03 0.045]);
% set_x0 = [X1(:) X2(:)];
% 
% set_x0 = [-2.5*ones(size(set_x0,1),1) set_x0 zeros(size(set_x0,1),1)];
% set_x0(:,end+1) = track.rref*ones(size(set_x0,1),1);
% 
% set_x0 = [0 0 0 0.025 track.rref];
% track.param_estimation = outputForParameterEstimation.empty;
% tic;
% for ct=1:1%size(set_x0,1)
%     
%     x0 = set_x0(ct,:)';
%     
%     [xOpt,cost,eflag,output] = patternsearch(@(x) landingDynamics_Objective(x(1), x(2), x(3), x(4), x(5), params),...
%         x0,[],[],[],[],params.lb([1 2 3 4 5]),params.ub([1 2 3 4 5]),[],options); % [0 0 0 0.025 0] [0 0 0 0.025 -1.8] [0.374 -7.4 -0.125 0.0263 -1.8]
%     
%     track.param_estimation(ct) = outputForParameterEstimation();
%     track.param_estimation(ct).x0 = x0;
%     track.param_estimation(ct).xOpt = xOpt;
%     track.param_estimation(ct).cost = cost;
%     track.param_estimation(ct).eflag = eflag;
%     track.param_estimation(ct).output = output;
% end
% toc;
% Kp = track.param_estimation(ct).xOpt(1); Ki = track.param_estimation(ct).xOpt(2); Kd = track.param_estimation(ct).xOpt(3); tau = track.param_estimation(ct).xOpt(4); rref = track.param_estimation(ct).xOpt(5);
% 
% figureHandles = [];

% Estimate Kp, Ki, Kd, tau
options = optimoptions('patternsearch','UseCompletePoll',true,'UseParallel',true,'Display','iter','MaxIterations',1000);

[X1, X2] = ndgrid([-5 0 5],[0.015 0.03 0.045]);
set_x0 = [X1(:) X2(:)];

set_x0 = [-2.5*ones(size(set_x0,1),1) set_x0 zeros(size(set_x0,1),1)];
set_x0(:,end+1) = track.rref*ones(size(set_x0,1),1);

set_x0 = [-5 5 0 0.025];
rref = track.rref;
track.param_estimation = outputForParameterEstimation.empty;
tic;
for ct=1:1%size(set_x0,1)
    
    x0 = set_x0(ct,:)';
    
    [xOpt,cost,eflag,output] = patternsearch(@(x) landingDynamics_Objective(x(1), x(2), x(3), x(4), rref, params),...
        x0,[],[],[],[],params.lb([1 2 3 4]),params.ub([1 2 3 4]),[],options); % [0 0 0 0.025 0] [0 0 0 0.025 -1.8] [0.374 -7.4 -0.125 0.0263 -1.8]
    
    track.param_estimation(ct) = outputForParameterEstimation();
    track.param_estimation(ct).x0 = x0;
    track.param_estimation(ct).xOpt = xOpt;
    track.param_estimation(ct).cost = cost;
    track.param_estimation(ct).eflag = eflag;
    track.param_estimation(ct).output = output;
end
toc;
Kp = track.param_estimation(ct).xOpt(1); Ki = track.param_estimation(ct).xOpt(2); Kd = track.param_estimation(ct).xOpt(3); tau = track.param_estimation(ct).xOpt(4);

figureHandles = [];


% Estimate Kp, Ki, Kd
options = optimoptions('patternsearch','UseCompletePoll',true,'UseParallel',true,'Display','iter','MaxIterations',1000);

[X1, X2] = ndgrid([-5 0 5],[0.015 0.03 0.045]);
set_x0 = [X1(:) X2(:)];

set_x0 = [-2.5*ones(size(set_x0,1),1) set_x0 zeros(size(set_x0,1),1)];
set_x0(:,end+1) = track.rref*ones(size(set_x0,1),1);

set_x0 = [-15 15 0];
tau = 0.025;
rref = track.rref;
track.param_estimation = outputForParameterEstimation.empty;
tic;
for ct=1:1%size(set_x0,1)
    
    x0 = set_x0(ct,:)';
    
    [xOpt,cost,eflag,output] = patternsearch(@(x) landingDynamics_Objective(x(1), x(2), x(3), tau, rref, params),...
        x0,[],[],[],[],params.lb([1 2 3]),params.ub([1 2 3]),[],options); % [0 0 0 0.025 0] [0 0 0 0.025 -1.8] [0.374 -7.4 -0.125 0.0263 -1.8]
    
    track.param_estimation(ct) = outputForParameterEstimation();
    track.param_estimation(ct).x0 = x0;
    track.param_estimation(ct).xOpt = xOpt;
    track.param_estimation(ct).cost = cost;
    track.param_estimation(ct).eflag = eflag;
    track.param_estimation(ct).output = output;
end
toc;
Kp = track.param_estimation(ct).xOpt(1); Ki = track.param_estimation(ct).xOpt(2); Kd = track.param_estimation(ct).xOpt(3);

figureHandles = [];


% Estimate Kp, Ki
options = optimoptions('patternsearch','UseCompletePoll',true,'UseParallel',true,'Display','iter','MaxIterations',1000);

[X1, X2] = ndgrid([-5 0 5],[0.015 0.03 0.045]);
set_x0 = [X1(:) X2(:)];

set_x0 = [-2.5*ones(size(set_x0,1),1) set_x0 zeros(size(set_x0,1),1)];
set_x0(:,end+1) = track.rref*ones(size(set_x0,1),1);

set_x0 = [5 -5];
Kd = 0;
tau = 0.030;
rref = track.rref;
track.param_estimation = outputForParameterEstimation.empty;
tic;
for ct=1:1%size(set_x0,1)
    
    x0 = set_x0(ct,:)';
    
    [xOpt,cost,eflag,output] = patternsearch(@(x) landingDynamics_Objective(x(1), x(2), Kd, tau, rref, params),...
        x0,[],[],[],[],params.lb([1 2]),params.ub([1 2]),[],options); % [0 0 0 0.025 0] [0 0 0 0.025 -1.8] [0.374 -7.4 -0.125 0.0263 -1.8]
    
    track.param_estimation(ct) = outputForParameterEstimation();
    track.param_estimation(ct).x0 = x0;
    track.param_estimation(ct).xOpt = xOpt;
    track.param_estimation(ct).cost = cost;
    track.param_estimation(ct).eflag = eflag;
    track.param_estimation(ct).output = output;
end
toc;
Kp = track.param_estimation(ct).xOpt(1); Ki = track.param_estimation(ct).xOpt(2); 

figureHandles = [];

%% Surrogateopt
% Estimating Kp, Ki, Kd, tau

[X1, X2, X3] = ndgrid([-5 0 5],[-5 0 5],[0.015 0.03 0.045]);
set_x0 = [X1(:) X2(:)];

set_x0 = [set_x0 zeros(size(set_x0,1),1) X3(:)];
rref = -3.335;

rng default
options = optimoptions('surrogateopt','MaxFunctionEvaluations',10000,'Display','iter',...
    'UseParallel',true,'MinSurrogatePoints',30,'InitialPoints',set_x0,'PlotFcn','surrogateoptplot');
problem = struct('objective',@(x) landingDynamics_Objective(x(1), x(2), x(3), x(4), rref, params),...
    'lb',params.lb([1 2 3 4]),...
    'ub',params.ub([1 2 3 4]),...
    'options',options,...
    'solver','surrogateopt');
[xOpt,cost,eflag,output,trials] = surrogateopt(problem);
Kp = xOpt(1); Ki = xOpt(2); Kd = xOpt(3); tau = xOpt(4);
%% Plots the best solution obtained

% params.Kp = x(1);
% params.Ki = x(2);
% params.Kd = x(3);
% params.tau = x(4);
% params.rref = x(5);

y_minustau =  history(-tau, params);
x = [Kp, Ki, Kd, tau, rref, y_minustau(1)]'; 
sol = ddensd(@(t,y,ydel,ypdel) ddefun(t,y,ydel,ypdel,x), ...
             tau, tau, @(t) history(t, params), params.tspan, params.ddeoptions);
[y,yp] = deval(sol,params.time); % to find state derivative at specified time instants

figure;
plot(y(1,:), y(2,:)); hold on;
plot(params.iodata(:,1), params.iodata(:,2),'r');
legend('Simulated','Measured');
xlabel('y_1');
ylabel('y_2');

figure;
subplot(3,1,1);
plot(params.time, y(1,:)); hold on;
plot(params.time, params.iodata(:,1), 'r');
ylabel('y_1');
legend('Simulated', 'Measured');
subplot(3,1,2);
plot(params.time, y(2,:)); hold on;
plot(params.time, params.iodata(:,2), 'r');
ylabel('y_2 or yp_1');
subplot(3,1,3);
plot(params.time, yp(2,:)); hold on;
plot(params.time, params.iodata(:,3), 'r');
xlabel('t');
ylabel('yp_2');

figure; % Plot prediction errors
residuals = ;


% 
% % Figure to check state and state derivatives in negative time
% negtime = [-tau:0.001:0]';
% y_negtime = history(negtime, params);
% figure;
% subplot(3,1,1);
% plot(negtime, y_negtime(1,:),'bo'); hold on;
% plot(params.negTime_state(:,1), params.negTime_state(:,2),'r');
% legend('Interpolated','Measured');
% subplot(3,1,2);
% plot(negtime, y_negtime(2,:),'bo'); hold on;
% plot(params.negTime_state(:,1), params.negTime_state(:,3),'r');
% subplot(3,1,3);
% plot(negtime, diffxy(negtime,y_negtime(2,:)),'bo'); hold on;
% plot(params.negTime_state(:,1), params.negTime_state(:,4),'r');




%% DEBUG - Check gradient evaluation
% x = [0.374 -7.4 -0.125 0.0263 -0.0312];
% Kp = x(1); Ki = x(2); Kd = x(3); tau = x(4); rref = x(5);
% 
% landingDynamics_Objective(Kp, Ki, Kd, tau, rref, params);
% compareGradientEvaluation(Kp, Ki, Kd, tau, rref, params);

% figureHandles(1) = figure;
% plot(time, iodata(:,1), 'r-'); hold on;
% plot(rSignal.Values.Time,rSignal.Values.Data,'o');
% legend('Measured r', 'Simulated r');
% legend(gca,'show','Location','best');
% ylabel('r (1/s)', 'FontSize', 16);
% xlabel('time (s)', 'FontSize', 16);
% set(gca, 'FontSize', 18); grid on;
% title(['Estimation using ' outputSignalForEstimation], 'FontSize', 12);
% 
% 
% figureHandles(2) = figure;
% p1 = subplot(2,1,1);
% plot(iodata(:,3), iodata(:,2), 'r-'); hold on;
% plot(ySignal.Values.Data,VgySignal.Values.Data,'s');
% legend('Measured', 'Simulated');
% legend(gca,'show','Location','best');
% ylabel('V_{gy} (m/s)', 'FontSize', 16);
% set(gca, 'FontSize', 18); grid on;
% title(['Estimation using ' outputSignalForEstimation], 'FontSize', 12);
% 
% 
% p2 = subplot(2,1,2);
% plot(iodata(:,3), iodata(:,4), 'r-'); hold on;
% plot(ySignal.Values.Data,aySignal.Values.Data,'^');
% % legend('Measured a_y', 'Simulated a_y');
% % legend(gca,'show','Location','best');
% ylabel('a_{y} (m/s^2)', 'FontSize', 16);
% xlabel('y (m)', 'FontSize', 16);
% set(gca, 'FontSize', 18); grid on;
% 
% linkaxes([p1 p2],'x');

% keyboard;

% plot(time, iodata, ...
%     rSignal.Values.Time,rSignal.Values.Data,'o', ...
%     VgySignal.Values.Time,VgySignal.Values.Data,'s', ...
%     ySignal.Values.Time,ySignal.Values.Data,'^');
% title('Simulated and Measured Responses After Estimation');
% legend('Measured r',  'Measured V_{gy}', 'Measured y', ...
%     'Simulated r', 'Simulated V_{gy}', 'Simulated y');
% 
% figureHandles(2) = figure;
% plot(time, iodata, ...
%     rSignal.Values.Time,rSignal.Values.Data,'o', ...
%     VgySignal.Values.Time,VgySignal.Values.Data,'s', ...
%     ySignal.Values.Time,ySignal.Values.Data,'^');
% title('Simulated and Measured Responses After Estimation');
% legend('Measured r',  'Measured V_{gy}', 'Measured y', ...
%     'Simulated r', 'Simulated V_{gy}', 'Simulated y');

end


function vals = landingDynamics_Objective(Kp, Ki, Kd, tau, rref, params)
%	 landingDynamicsObjective
%
%    The landingDynamics_Objective function is used to calculate the objective function and compare model
%    outputs against experimental data.
%
%    The |vals| return argument contains information about how well the
%    model simulation results match the experimental data and is used by
%    the optimizer to estimate the model parameters.
%

% Create params structure
% This param structure is fixed to evaluate a cost function for a fixed set
% of parameters. y_minustau needs to be updated with new set of parameters
% during optimization.

% For each call to objective function, construct x
y_minustau =  history(-tau, params); % y(-tau)
x = [Kp, Ki, Kd, tau, rref, y_minustau(1)]'; % parameters

% sol = ddensd(@(t,y,ydel,ypdel) ddefun(t,y,ydel,ypdel,x), ...
%              tau, tau, @(t) history(t, params), [params.tspan(1)-tau params.tspan(2)], params.ddeoptions);

sol = ddensd(@(t,y,ydel,ypdel) ddefun(t,y,ydel,ypdel,x), ...
             tau, tau, @(t) history(t, params), params.tspan, params.ddeoptions);
         
% sol = ddensd(@(t,y,ydel,ypdel) ddefun(t,y,ydel,ypdel,x), ...
%              tau, tau, params.histsol, [params.tspan(1)-tau params.tspan(2)], params.ddeoptions);
% 
% sol = ddensd(@(t,y,ydel,ypdel) ddefun(t,y,ydel,ypdel,x), ...
%              tau, tau, params.histsol, params.tspan, params.ddeoptions);
         
if params.tspan(2)>sol.x(end) % Event specified in ddeoptions has occured 
    vals = 200; % Return some arbitrarily high value of the cost function
    return
%     keyboard;
end
         
[y,yp] = deval(sol,params.time); % to find state derivative at specified time instants

% Evaluating cost function
if strcmpi(params.cost, 'sae') % and algorithm is not lsqnonlin
    % Weighted SAE
%     vals = sum(params.weights.*abs(yp(2,:)'-params.iodata(:,3))); % accleration
%     vals = sum(params.weights.*abs(y(1,:)'-params.iodata(:,1))); % position
    vals = sum(params.weights.*abs(y(2,:)'./y(1,:)'-params.iodata(:,2)./params.iodata(:,1))); % vel/position
elseif strcmpi(params.cost, 'sse')
    % Weighted SSE
    vals = sum(params.weights.*(yp(2,:)'-params.iodata(:,3)).^2);
else
    error('Undefined cost function..');
end

% figure;
% plot(y(1,:), y(2,:)); hold on;
% plot(params.iodata(:,1), params.iodata(:,2),'r');
% legend('Simulated','Measured');
% xlabel('y_1');
% ylabel('y_2');
% 
% figure;
% plot(params.time, yp(2,:)); hold on;
% plot(params.time, params.iodata(:,3), 'r');
% legend('Simulated', 'Measured');
% xlabel('t');
% ylabel('yp_2');

end

function y = history(t,params)
% Provides y vector for negative times
% y is a state vector (=[y; Vgy])
% params - a custom structure containing data that remain fixed within one call to estimateLandingParameters.m

y = interp1(params.negTime_state(:,1), params.negTime_state(:,2:3), t, 'linear')';
% y = params.y0;
end

function yp = ddefun(t,y,ydel,ypdel,x)
% Set of neutral-delay-differential equations governing closed-loop landing
% dynamics of bbees

% x is a vector of parameters that can change for different calls to the
% objective function evaluation (= [Kp Ki Kd tau rref y_minustau]), where
% y_minustau is a y-coordinate (first element of state vector y) at t=-tau

% y is a state vector (=[y; Vgy])
% yp is dy/dt
% ydel and ypdel contain delays in state and state derivatives
% ydel and ypdel have only one column corresponding to tau (sensory delay)

y1dot = y(2);
y2dot = x(1)*(x(5)-ydel(2)/ydel(1)) + ...
        x(2)*(x(5)*t-log(abs(ydel(1)))+log(abs(x(6)))) - ...
        x(3)*(ypdel(2)./ydel(1)-(ydel(2)./ydel(1)).^2);

yp = [y1dot;
      y2dot];
end

function [position,isterminal,direction] = events(t,y,ydel,ypdel)
% Terminate when bbee hits the platform
position = y(1); % Detect y = 0
isterminal = 1; % Stop the integration
direction = 1; % Positive direction only (because y=0 is defined at the platform and y-axis goes inside the platform
end

function [cf, cfgrad] = landingDynamics_scalarObjective_Jacobian(Kp, Ki, Kd, tau, rref, params)
%	 landingDynamics_scalarObjective_Jacobian
%
%    The |vals| return argument contains information about how well the
%    model simulation results match the experimental data and is used by
%    the optimizer to estimate the model parameters.
%

% Create params structure
% This param structure is fixed to evaluate a cost function for a fixed set
% of parameters. y_minustau needs to be updated with new set of parameters
% during optimization.

% For each call to objective function, construct x
y_minustau =  history(-tau, params); % y(-tau)
x = [Kp, Ki, Kd, tau, rref, y_minustau(1)]'; % parameters

sol = ddensd(@(t,y,ydel,ypdel) ddefun(t,y,ydel,ypdel,x), ...
             tau, tau, @(t) history(t, params), [params.tspan(1)-tau-params.epsilon params.tspan(2)], params.ddeoptions);

if params.tspan(2)>sol.x(end) % Event specified in ddeoptions has occured 
    cf = 1000; % Return some arbitrarily high value of the cost function
    cfgrad = [Inf Inf Inf Inf Inf];
    return
%     keyboard;
end
 
[~,yp_t] = deval(sol,params.time); % to find state derivative at specified time instants
[y_tminustau,yp_tminustau] = deval(sol,params.time-tau); % to find state derivative at specified time instants

[~,yp_t_minustau_minusepsilon] = deval(sol,params.time-tau-params.epsilon); % to find state derivative at specified time instants
[~,yp_t_minustau_plusepsilon] = deval(sol,params.time-tau+params.epsilon); % to find state derivative at specified time instants
ypp_tminustau = (yp_t_minustau_plusepsilon-yp_t_minustau_minusepsilon)./(2*params.epsilon);
% Evaluating cost function
if strcmpi(params.cost, 'sae')
    % Weighted SAE
    cf = sum(params.weights.*abs(yp_t(2,:)'-params.iodata(:,3)));
    dJ_dKp = sum(-1*params.weights.*sign(params.iodata(:,3)-yp_t(2,:)').*(rref-y_tminustau(2,:)'./y_tminustau(1,:)'));
    dJ_dKi = sum(-1*params.weights.*sign(params.iodata(:,3)-yp_t(2,:)').*(rref*params.time-log(abs(y_tminustau(1,:)'))+log(abs(y_minustau(1)))));
    dJ_dKd = sum(params.weights.*sign(params.iodata(:,3)-yp_t(2,:)').*(yp_tminustau(2,:)'./y_tminustau(1,:)'-(y_tminustau(2,:)./y_tminustau(1,:))'.^2));
    dJ_dtau = sum(-1*params.weights.*sign(params.iodata(:,3)-yp_t(2,:)').*( ...
        Kp*(-(y_tminustau(2,:)./y_tminustau(1,:))'.^2 + yp_tminustau(2,:)'./y_tminustau(1,:)') ...
       +Ki*(y_tminustau(2,:)'./y_tminustau(1,:)' - y_minustau(2)/y_minustau(1)) ...
       +Kd*(-(yp_tminustau(2,:)./y_tminustau(1,:))'.^2 + ypp_tminustau(2,:)'./y_tminustau(1,:)' ...
            + 2*(y_tminustau(2,:)./y_tminustau(1,:))'.^3 - 2*yp_tminustau(2,:)'.*y_tminustau(2,:)'./y_tminustau(1,:)'.^2) ...
       ));
    dJ_drref = sum(-1*params.weights.*sign(params.iodata(:,3)-yp_t(2,:)').*(Kp+Ki*params.time));
    cfgrad = [dJ_dKp dJ_dKi dJ_dKd dJ_dtau dJ_drref];
elseif strcmpi(params.cost, 'sse')
    % Weighted SSE
    cf = sum(params.weights.*(yp_t(2,:)'-params.iodata(:,3)).^2);
    cfgrad = [];
else
    error('Undefined cost function..');
end

end

function compareGradientEvaluation(Kp, Ki, Kd, tau, rref, params)
% Analytical computation
[cf, cfgrad] = landingDynamics_scalarObjective_Jacobian(Kp, Ki, Kd, tau, rref, params);

disp('Cost functions:');
disp([cf landingDynamics_Objective(Kp, Ki, Kd, tau, rref, params)]);

disp('Analytically computed gradient vector:');
disp(num2str(cfgrad));

epsilon = 1e-2;
% forward difference estimate
cfgrad_fd = [ (landingDynamics_Objective(Kp+epsilon, Ki, Kd, tau, rref, params)-landingDynamics_Objective(Kp, Ki, Kd, tau, rref, params))/epsilon, ...
              (landingDynamics_Objective(Kp, Ki+epsilon, Kd, tau, rref, params)-landingDynamics_Objective(Kp, Ki, Kd, tau, rref, params))/epsilon, ...
              (landingDynamics_Objective(Kp, Ki, Kd+epsilon, tau, rref, params)-landingDynamics_Objective(Kp, Ki, Kd, tau, rref, params))/epsilon, ...
              (landingDynamics_Objective(Kp, Ki, Kd, tau+epsilon, rref, params)-landingDynamics_Objective(Kp, Ki, Kd, tau, rref, params))/epsilon, ...
              (landingDynamics_Objective(Kp, Ki, Kd, tau, rref+epsilon, params)-landingDynamics_Objective(Kp, Ki, Kd, tau, rref, params))/epsilon];

disp('Forward difference estimate of gradient vector:');
disp(num2str(cfgrad_fd));
% Central difference estimate
cfgrad_cd = [ (landingDynamics_Objective(Kp+epsilon, Ki, Kd, tau, rref, params)-landingDynamics_Objective(Kp-epsilon, Ki, Kd, tau, rref, params))/epsilon/2, ...
              (landingDynamics_Objective(Kp, Ki+epsilon, Kd, tau, rref, params)-landingDynamics_Objective(Kp, Ki-epsilon, Kd, tau, rref, params))/epsilon/2, ...
              (landingDynamics_Objective(Kp, Ki, Kd+epsilon, tau, rref, params)-landingDynamics_Objective(Kp, Ki, Kd-epsilon, tau, rref, params))/epsilon/2, ...
              (landingDynamics_Objective(Kp, Ki, Kd, tau+epsilon, rref, params)-landingDynamics_Objective(Kp, Ki, Kd, tau-epsilon, rref, params))/epsilon/2, ...
              (landingDynamics_Objective(Kp, Ki, Kd, tau, rref+epsilon, params)-landingDynamics_Objective(Kp, Ki, Kd, tau, rref-epsilon, params))/epsilon/2];

disp('Central difference estimate of gradient vector:');
disp(num2str(cfgrad_cd));
end

function [cf, cfgrad] = DEBUG_landingDynamics_scalarObjective_Jacobian(Kp, Ki, Kd, tau, rref, params)
%	 landingDynamics_scalarObjective_Jacobian
%
%    The |vals| return argument contains information about how well the
%    model simulation results match the experimental data and is used by
%    the optimizer to estimate the model parameters.
%

% Create params structure
% This param structure is fixed to evaluate a cost function for a fixed set
% of parameters. y_minustau needs to be updated with new set of parameters
% during optimization.

% For each call to objective function, construct x
y_minustau =  history(-tau, params); % y(-tau)
x = [Kp, Ki, Kd, tau, rref, y_minustau(1)]'; % parameters

sol = ddensd(@(t,y,ydel,ypdel) ddefun(t,y,ydel,ypdel,x), ...
             tau, tau, @(t) history(t, params), [params.tspan(1)-tau-params.epsilon params.tspan(2)], params.ddeoptions);

if params.tspan(2)>sol.x(end) % Event specified in ddeoptions has occured 
    cf = 1000; % Return some arbitrarily high value of the cost function
    cfgrad = [Inf Inf Inf Inf Inf];
    return
%     keyboard;
end
 
[~,yp_t] = deval(sol,params.time); % to find state derivative at specified time instants
[y_tminustau,yp_tminustau] = deval(sol,params.time-tau); % to find state derivative at specified time instants

[~,yp_t_minustau_minusepsilon] = deval(sol,params.time-tau-params.epsilon); % to find state derivative at specified time instants
[~,yp_t_minustau_plusepsilon] = deval(sol,params.time-tau+params.epsilon); % to find state derivative at specified time instants
ypp_tminustau = (yp_t_minustau_plusepsilon-yp_t_minustau_minusepsilon)./(2*params.epsilon);
% Evaluating cost function
if strcmpi(params.cost, 'sae')
    % Weighted SAE
    cf = sum(params.weights.*abs(yp_t(2,:)'-params.iodata(:,3)));
    dJ_dKp = sum(-1*params.weights.*sign(params.iodata(:,3)-yp_t(2,:)').*(rref-y_tminustau(2,:)'./y_tminustau(1,:)'));
    dJ_dKi = sum(-1*params.weights.*sign(params.iodata(:,3)-yp_t(2,:)').*(rref*params.time-log(abs(y_tminustau(1,:)'))+log(abs(y_minustau(1)))));
    dJ_dKd = sum(params.weights.*sign(params.iodata(:,3)-yp_t(2,:)').*(yp_tminustau(2,:)'./y_tminustau(1,:)'-(y_tminustau(2,:)./y_tminustau(1,:))'.^2));
    dJ_dtau = sum(-1*params.weights.*sign(params.iodata(:,3)-yp_t(2,:)').*( ...
        Kp*(-(y_tminustau(2,:)./y_tminustau(1,:))'.^2 + yp_tminustau(2,:)'./y_tminustau(1,:)') ...
       +Ki*(y_tminustau(2,:)'./y_tminustau(1,:)' - y_minustau(2)/y_minustau(1)) ...
       +Kd*(-(yp_tminustau(2,:)./y_tminustau(1,:))'.^2 + ypp_tminustau(2,:)'./y_tminustau(1,:)' ...
            + 2*(y_tminustau(2,:)./y_tminustau(1,:))'.^3 - 2*yp_tminustau(2,:)'.*y_tminustau(2,:)'./y_tminustau(1,:)'.^2) ...
       ));
    dJ_drref = sum(-1*params.weights.*sign(params.iodata(:,3)-yp_t(2,:)').*(Kp+Ki*params.time));
    cfgrad = [dJ_dKp dJ_dKi dJ_dKd];% dJ_dtau dJ_drref];
elseif strcmpi(params.cost, 'sse')
    % Weighted SSE
    cf = sum(params.weights.*(yp_t(2,:)'-params.iodata(:,3)).^2);
    cfgrad = [];
else
    error('Undefined cost function..');
end

end