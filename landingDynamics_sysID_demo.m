function landingDynamic_sysID_demo()
% LANDINGPARAMETERESTIMTION
%   Learns parameters of the landing model based on the landingModel.slx
%   track - instance of trackForLandingModel
%   runDemo - bool (if demo has to run)

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

[time,~,iodata] = sim('landingDynamics', 2);

%%
% Create an experiment object to store the measured data
Exp = sdo.Experiment('landingDynamics');

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

% Defining parameters to be estimated and setting values of some parameters
p = sdo.getParameterFromModel('landingDynamics',{'r_ref','output_delay','Kp', 'Ki', 'Kd', 'v0', 'y0'});

% setting scale of all parameters
[p.Scale] = deal(1);

% r_ref
p(1).Value = 0;
p(1).Minimum = -10;
p(1).Maximum = 0;
p(1).Free = true;

% output_delay
p(2).Value = 0.01;
p(2).Minimum = 0.001;
p(2).Maximum = 0.1;
p(2).Free = true;

% Kp
p(3).Value = 0;
p(3).Free = true;

% Ki
p(4).Value = 0;
p(4).Free = true;

% Kd
p(5).Value = 0;
p(5).Free = true;

% v0
p(6).Value = iodata(1,2);
p(6).Free = false;

% y0
p(7).Value = iodata(1,3);
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

%  defining parameters to be estimated
v = p(1:3)';

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
opt.MethodOptions = optimoptions('patternsearch','SearchFcn',{@searchga, 3}, ...
    'PlotFcn','psplotbestf', 'FunctionTolerance', 0.001, 'Display', 'iter');%, ...
%     'UseParallel', true, 'UseCompletePoll', true, 'UseCompleteSearch', true);

% Defining objective function
estFcn = @(v) landingDynamics_Objective(v, Simulator, Exp);

% Estimating the parameters
[vOpt, opt_info] = sdo.optimize(estFcn,v,opt);

%%
% Update the experiment with estimated parameter values
Exp = setEstimatedValues(Exp,vOpt);

% Compare measured and simulated data
Simulator    = createSimulator(Exp,Simulator);
Simulator    = sim(Simulator);

SimLog       = find(Simulator.LoggedData,get_param('landingDynamics','SignalLoggingName'));
rSignal      = find(SimLog,'r');
VgySignal    = find(SimLog,'Vgy');
ySignal      = find(SimLog,'y');

plot(time, iodata, ...
    rSignal.Values.Time,rSignal.Values.Data,'-.', ...
    VgySignal.Values.Time,VgySignal.Values.Data,'--', ...
    ySignal.Values.Time,ySignal.Values.Data,'-*');
title('Simulated and Measured Responses After Estimation')
legend('Measured r',  'Measured V_{gy}', 'Measured y', ...
    'Simulated r', 'Simulated V_{gy}', 'Simulated y');

end

function vals = landingDynamics_Objective(v,Simulator,Exp) 
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
ref_sig = sdo.requirements.SignalTracking;
ref_sig.Type      = '==';
% ref_sig.Method    = 'Residuals'; % can also be SSE or SAE
ref_sig.Method    = 'SAE'; % can also be SSE or SAE
ref_sig.Normalize = 'off';

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

rError       = evalRequirement(ref_sig, rSignal.Values, Exp.OutputData(1).Values);

%%
% Return the residual errors to the optimization solver.
%
vals.F = [rError(:)];
end
