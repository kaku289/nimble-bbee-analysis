%% This is an optional file to run the parametric sweep in order to evaluate the cost function.
% The cost function evaluation in this script is based on the simulink
% model and needs to be updated with the system of
% neutral-delayed-differential equations. Therefore, it is not recommended
% to run this program in this form. It took more than 30 hours last time to
% make ~33500 functional evaluations.
%%
clc; clear;

% 1) Load model
model = 'landingDynamics';
load_system(model);

% 2) Set up the sweep parameters
% Kp_sweep = -5:0.25:5;
% Ki_sweep = -5:0.25:5;
% Kd_sweep = -1:0.2:1;
% tau_sweep = 0.01:0.002:0.06;
% rref_sweep = -2:0.25:2;

Kp_sweep = -5:1:5;
Ki_sweep = -5:1:5;
Kd_sweep = -1:0.5:1;
tau_sweep = 0.01:0.005:0.06;
rref_sweep = -2:0.5:0;

numSims   = numel(Kp_sweep)*numel(Ki_sweep)*numel(Kd_sweep)*numel(tau_sweep)*numel(rref_sweep);

% 3) Create an array of SimulationInput objects and specify the sweep value for each simulation
simIn(1:numSims) = Simulink.SimulationInput(model);
% [x1, x2, x3, x4, x5] = ndgrid(Kp_sweep,Ki_sweep,Kd_sweep,tau_sweep,rref_sweep);

ct=1;
for ct1 = 1:numel(Kp_sweep)
    for ct2 = 1:numel(Ki_sweep)
        for ct3 = 1:numel(Kd_sweep)
            for ct4 = 1:numel(tau_sweep)
                for ct5 = 1:numel(rref_sweep)
                    simIn(ct) = simIn(ct).setVariable('Kp',Kp_sweep(ct1),'Workspace','landingDynamics');
                    simIn(ct) = simIn(ct).setVariable('Ki',Ki_sweep(ct2),'Workspace','landingDynamics');
                    simIn(ct) = simIn(ct).setVariable('Kd',Kd_sweep(ct3),'Workspace','landingDynamics');
                    simIn(ct) = simIn(ct).setVariable('output_delay',tau_sweep(ct4),'Workspace','landingDynamics');
                    simIn(ct) = simIn(ct).setVariable('r_ref',rref_sweep(ct5),'Workspace','landingDynamics');
                                                      
                    ct=ct+1;
                end
            end
        end
    end
end

% 4) Simulate the model 
simOut = parsim(simIn);

% Loading the reference signal
load('/home/reken001/Pulkit/MATLAB/graphs/checkerboard_high_t1e1_refSignal.mat');

% assigning weights
weights = ones(length(ref_signal.Time),1);

% Initializing variables
params = zeros(length(simOut),5);
cf_sae = zeros(length(params),1);
cf_sse = cf_sae;

% Define a signal tracking requirement
ref_sig_sae = sdo.requirements.SignalTracking('ReferenceSignal', ref_signal ,'Weights', weights);
ref_sig_sae.Type      = '==';
ref_sig_sae.Method    = 'SAE'; % can also be SSE or SAE
% ref_sig.Normalize = 'off';

ref_sig_sse = sdo.requirements.SignalTracking('ReferenceSignal', ref_signal ,'Weights', weights);
ref_sig_sse.Type      = '==';
ref_sig_sse.Method    = 'SSE'; % can also be SSE or SAE
% ref_sig.Normalize = 'off';

% Compute cost function
ct = 1;
for ct1 = 1:numel(Kp_sweep)
    for ct2 = 1:numel(Ki_sweep)
        for ct3 = 1:numel(Kd_sweep)
            for ct4 = 1:numel(tau_sweep)
                for ct5 = 1:numel(rref_sweep)
                    
                    % Saving parameter set
                    params(ct,:) = [Kp_sweep(ct1), Ki_sweep(ct2), Kd_sweep(ct3), tau_sweep(ct4), rref_sweep(ct5)];
                    
                    % cost functions
                    rSignal      = simOut(ct).yout{4};
                    cf_sae(ct)   = evalRequirement(ref_sig_sae, rSignal.Values);
                    cf_sse(ct)   = evalRequirement(ref_sig_sse, rSignal.Values);
                    
                    %
                    ct = ct + 1;
                end
            end
        end
    end
end

% Minimum cost function
disp('SAE cost function')
[min_sae, indx_sae] = min(cf_sae);
disp(['Minimum cf_sae: ' num2str(min_sae) ', Parameters: ' num2str(params(indx_sae,:))]);

disp('SSE cost function')
[min_sse, indx_sse] = min(cf_sse);
disp(['Minimum cf_sse: ' num2str(min_sse) ', Parameters: ' num2str(params(indx_sse,:))]);
