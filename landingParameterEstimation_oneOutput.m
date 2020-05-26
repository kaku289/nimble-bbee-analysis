function [vOpt, opt_info, params, figureHandles] = landingParameterEstimation_oneOutput(track, runDemo, outputSignalForEstimation)

%   LANDINGPARAMETERESTIMTIONONEOUTPUT
%   Learns parameters of the landing model based on the landingModel.slx
%   track - instance of trackForLandingModel
%   runDemo - bool (if demo has to run)
%   outputSignalForEstimation - 'r' or 'ay'


% SDOAIRCRAFTESTIMATION_EXPERIMENT
%
%    [time,iodata] = sdoAircraftEstimation_Experiment 
%
%    The sdoAircraftEstimation_Experiment function is used create experiment
%    data from the sdoAircraftEstimation model.
%
%    The |time| return argument is a vector of time points at which
%    experimental data is collected. The |iodata| return argument contains
%    experiment output data; the 1st column is the angle-of-attack signal
%    and the 2nd column the G-force experienced by the pilot.
%
%    See also sdoAircraftExperiment_cmddemo
%

%%
open_system('landingDynamics');

% Put measured data in form suitable for system identification
if runDemo
    mdlWks = get_param('landingDynamics','ModelWorkspace');
    
    % setting values of v0 and y0
    assignin(mdlWks, 'v0', 0.7);
    assignin(mdlWks, 'y0', -1);
    
    % setting values of other parameters
    assignin(mdlWks, 'r_ref', -0.5);
    assignin(mdlWks, 'output_delay', 0.05);
    assignin(mdlWks, 'Kp', -20);
    assignin(mdlWks, 'Ki', 0);
    assignin(mdlWks, 'Kd', 0);
    
    % setting value of negative time vs r_true
    assignin(mdlWks, 'negTime_r', [0 -0.5; -1 -0.5]);
    
    [time,~,iodata] = sim('landingDynamics', 2);
    negTime_r = [0 -0.5; -1 -0.5];
else
    time = track.state(:,1)-track.state(1,1);
    iodata = [track.state(:,6)./track.state(:,3), track.state(:,6), track.state(:,3), track.state(:,9)]; % [r Vgy y ay]
    negTime_r = [track.negTime_r(:,1)-track.negTime_r(end,1) track.negTime_r(:,6)./track.negTime_r(:,3)];
    v0 = iodata(1,2);
    y0 = iodata(1,3);
end


% Constructing weights signal
weights = ones(length(time),1);
% weights(time>(time(end)-0.1)) = 2;

%%
% Create an experiment object to store the measured data
Exp = sdo.Experiment('landingDynamics');

if strcmpi(outputSignalForEstimation, 'r')
    % Creating r (true optical rate of expansion - as obtained from the
    % processing of the trajectories)
    r = Simulink.SimulationData.Signal;
    r.Name = 'r';
    r.BlockPath = 'landingDynamics/Open-loop dynamics';
    r.PortType = 'outport';
    r.PortIndex = 1;
    r.Values = timeseries(iodata(:,1), time);

    % Adding measured data to the experiment
    Exp.OutputData = [r];
    
elseif strcmpi(outputSignalForEstimation, 'ay')
    % Creating r (true optical rate of expansion - as obtained from the
    % processing of the trajectories)
    ay = Simulink.SimulationData.Signal;
    ay.Name = 'r';
    ay.BlockPath = 'landingDynamics/Open-loop dynamics';
    ay.PortType = 'outport';
    ay.PortIndex = 4;
    ay.Values = timeseries(iodata(:,4), time);

    % Adding measured data to the experiment
    Exp.OutputData = [ay];
else
    error('Unknown string for outputSignalForEstimation is found. It can be either r or ay.');
end


% Defining parameters to be estimated and setting values of some parameters
p = sdo.getParameterFromModel('landingDynamics',{'r_ref','output_delay','Kp', 'Ki', 'Kd', 'v0', 'y0'});

% setting scale of all parameters
% [p.Scale] = deal(1);

% r_ref
p(1).Value = 0;
p(1).Minimum = -10;
p(1).Maximum = 10;
p(1).Free = true;

% output_delay
p(2).Value = 0.020;
p(2).Minimum = 0.01;
p(2).Maximum = 0.06;
% p(2).Scale = 0.0156;
p(2).Free = true;

% Kp
p(3).Value = 0;
p(3).Free = true;
p(3).Minimum = -10;
p(3).Maximum = 10;

% Ki
p(4).Value = 0;
p(4).Free = true;
% p(4).Free = false;
p(4).Minimum = -10;
p(4).Maximum = 10;

% Kd
p(5).Value = 0;
p(5).Free = true;
% p(5).Free = false;
p(5).Minimum = -10;
p(5).Maximum = 10;

% v0
p(6).Value = v0; %iodata(1,2);
p(6).Free = false;

% y0
p(7).Value = y0; %iodata(1,3);
p(7).Free = false;

% Get model workspace of landingDynamics.slx
mdlWks = get_param('landingDynamics','ModelWorkspace');

% setting values of v0 and y0
assignin(mdlWks, 'v0', p(6).Value);
assignin(mdlWks, 'y0', p(7).Value);

% setting values of each parameter
assignin(mdlWks, 'r_ref', p(1).Value);
assignin(mdlWks, 'output_delay', p(2).Value);
assignin(mdlWks, 'Kp', p(3).Value);
assignin(mdlWks, 'Ki', p(4).Value);
assignin(mdlWks, 'Kd', p(5).Value);

% setting value of negative time vs r_true
assignin(mdlWks, 'negTime_r', negTime_r);

%  defining parameters to be estimated
if runDemo
    v = p(1:3);
else
    v = p(1:5);
end

%%
% Creating a simulator
Simulator    = createSimulator(Exp);
Simulator    = sim(Simulator);

%% 
% Defining optimization options
% https://nl.mathworks.com/help/sldo/ref/sdo.optimizeoptions-class.html
% https://nl.mathworks.com/help/optim/ug/optim.problemdef.optimizationproblem.optimoptions.html
% https://nl.mathworks.com/help/gads/pattern-search-options.html

opt = sdo.OptimizeOptions;

% lsqnonlin method
% opt.Method = 'lsqnonlin';

% pattern search method
opt.Method = 'patternsearch';
a = load('/home/reken001/Pulkit/MATLAB/landingDynamics/data4estimationUsingGUI_spesession.mat');
opt.MethodOptions = a.SDOSessionData.Data.ToolData.ParameterEstimation.Options.MethodOptions;
opt.MethodOptions.MaxIter = 200;
opt.MethodOptions.CompletePoll = 'on';
% opt.MethodOptions.UseParallel = true;
% opt.MethodOptions = optimoptions('patternsearch','SearchFcn','searchga', ...
%     'PlotFcn','psplotbestf', 'FunctionTolerance', 0.001, 'Display', 'iter', 'MaxIterations', 2*length(v)*100, ...
%     'UseCompleteSearch', true);
% opt.MethodOptions = optimoptions('patternsearch','SearchFcn',{@searchga, 3}, ...
%     'PlotFcn','psplotbestf', 'FunctionTolerance', 0.001, 'Display', 'iter', 'MaxIterations', length(v)*100);%, ...
%     'UseParallel', true, 'UseCompletePoll', true, 'UseCompleteSearch', true);

% Defining objective function
estFcn = @(v) landingDynamics_Objective(v, Simulator, Exp, weights);

% Estimating the parameters
[vOpt, opt_info] = sdo.optimize(estFcn,v,opt);

%%
% Update the experiment with estimated parameter values
Exp = setEstimatedValues(Exp,vOpt);

% Making updated set of parameter values (for output)
params = [vOpt; p(length(vOpt)+1:end)];

% Creating r 
r = Simulink.SimulationData.Signal;
r.Name = 'r';
r.BlockPath = 'landingDynamics/Open-loop dynamics';
r.PortType = 'outport';
r.PortIndex = 1;
r.Values = timeseries(iodata(:,1), time);

% Creating Vgy
Vgy = Simulink.SimulationData.Signal;
Vgy.Name = 'Vgy';
Vgy.BlockPath = 'landingDynamics/Open-loop dynamics';
Vgy.PortType = 'outport';
Vgy.PortIndex = 2;
Vgy.Values = timeseries(iodata(:,2), time);

% Creating y
y = Simulink.SimulationData.Signal;
y.Name = 'y';
y.BlockPath = 'landingDynamics/Open-loop dynamics';
y.PortType = 'outport';
y.PortIndex = 3;
y.Values = timeseries(iodata(:,3), time);

% Creating ay
ay = Simulink.SimulationData.Signal;
ay.Name = 'ay';
ay.BlockPath = 'landingDynamics/Open-loop dynamics';
ay.PortType = 'outport';
ay.PortIndex = 4;
ay.Values = timeseries(iodata(:,4), time);

% Define output data of the experiment
Exp.OutputData = [r; Vgy; y; ay];

% Compare measured and simulated data
Simulator    = createSimulator(Exp);
% Simulator    = createSimulator(Exp,Simulator);
Simulator    = sim(Simulator);

SimLog       = find(Simulator.LoggedData,get_param('landingDynamics','SignalLoggingName'));
rSignal      = find(SimLog,'r');
VgySignal    = find(SimLog,'Vgy');
ySignal      = find(SimLog,'y');
aySignal     = find(SimLog,'ay');

figureHandles(1) = figure;
plot(time, iodata(:,1), 'r-'); hold on;
plot(rSignal.Values.Time,rSignal.Values.Data,'o');
legend('Measured r', 'Simulated r');
legend(gca,'show','Location','best');
ylabel('r (1/s)', 'FontSize', 16);
xlabel('time (s)', 'FontSize', 16);
set(gca, 'FontSize', 18); grid on;
title(['Estimation using ' outputSignalForEstimation], 'FontSize', 12);


figureHandles(2) = figure;
p1 = subplot(2,1,1);
plot(iodata(:,3), iodata(:,2), 'r-'); hold on;
plot(ySignal.Values.Data,VgySignal.Values.Data,'s');
legend('Measured', 'Simulated');
legend(gca,'show','Location','best');
ylabel('V_{gy} (m/s)', 'FontSize', 16);
set(gca, 'FontSize', 18); grid on;
title(['Estimation using ' outputSignalForEstimation], 'FontSize', 12);


p2 = subplot(2,1,2);
plot(iodata(:,3), iodata(:,4), 'r-'); hold on;
plot(ySignal.Values.Data,aySignal.Values.Data,'^');
% legend('Measured a_y', 'Simulated a_y');
% legend(gca,'show','Location','best');
ylabel('a_{y} (m/s^2)', 'FontSize', 16);
xlabel('y (m)', 'FontSize', 16);
set(gca, 'FontSize', 18); grid on;

linkaxes([p1 p2],'x');

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


function vals = landingDynamics_Objective(v,Simulator,Exp, weights) 
%	 landingDynamicsObjective
%
%    The landingDynamics_Objective function is used to calculate the objective function and compare model
%    outputs against experimental data.
%
%    vals = landingDynamics_Objective(v,Simulator,Exp)
%
%    The |v| input argument is a vector of estimated model parameter values
%    (and initial states, if any).
%
%    The |Simulator| input argument is a simulation object used 
%    simulate the model with the estimated parameter values.
%
%    The |Exp| input argument contains the estimation experiment data.
%
%    The |vals| return argument contains information about how well the
%    model simulation results match the experimental data and is used by
%    the |sdo.optimize| function to estimate the model parameters.
%

%%
% Define a signal tracking requirement to compute how well the model output
% matches the experiment data. Configure the tracking requirement so that
% it returns the tracking error residuals (rather than the
% sum-squared-error, or sum-absolute errors) and does not normalize the errors.
% https://nl.mathworks.com/help/sldo/ref/sdo.requirements.signaltracking-class.html
ref_sig = sdo.requirements.SignalTracking('ReferenceSignal', Exp.OutputData(1).Values ,'Weights', weights);
ref_sig.Type      = '==';
% ref_sig.Method    = 'Residuals'; % can also be SSE or SAE
ref_sig.Method    = 'SAE'; % can also be SSE or SAE
% ref_sig.Normalize = 'off';

%%
% Update the experiments with the estimated parameter values.
% https://nl.mathworks.com/help/sldo/ref/sdo.experiment.setestimatedvalues.html
Exp  = setEstimatedValues(Exp,v);

%%
% Simulate the model and compare model outputs with measured experiment
% data.
%
Simulator = createSimulator(Exp,Simulator);
Simulator = sim(Simulator);

SimLog       = find(Simulator.LoggedData,get_param('landingDynamics','SignalLoggingName'));
rSignal      = find(SimLog,'r');

% rError       = evalRequirement(ref_sig, rSignal.Values, Exp.OutputData(1).Values);
rError       = evalRequirement(ref_sig, rSignal.Values);

%%
% Return the errors to the optimization solver.
%
vals.F = [rError(:)];
end
