classdef trackForLandingModel < handle
    properties (Access = public)
        y_start = 0;
        y_end = 0;
        
        
        state = []; % storing state to be used for optimization (for positive times)
        negTime_state = []; % To store state for negative time
        
        % values set after optimization is complete
        param_estimation = outputForParameterEstimation.empty; % N by 1 array, where N is the number of initial conditions
        
        % state for rref estimation
        data4rrefEstimate = data4rrefEstimate.empty; % n_part by 1 array, where n_part is the number of identified segments within each track excerpt with different rref
        
%         state4rrefEstimate = []; % NX10 vector - the whole shebang
%         rref = nan;
%         vmean = nan; % Mean of Vgy in state4rrefEstimate
%         ymean = nan; % Mean of y in state4rrefEstimate
%         model; % linear fit model
    end
     methods
         function obj = trackForLandingModel()
             
         end
         
         function plotHandles = estimate_rref(object, createPlots, state_LDF)
             % state_LDF - The whole state from which state4rrefEstimate is
             % extracted (given as an instance of filteredState_BlindLandingtrack class)
             
             % plotHandles(1) - Plotted with state_LDF
             % plotHandles(2) - Only the data that is used to compute rref (code needs to be altered to produce it correctly)
             if ~isempty(object.data4rrefEstimate)
                 
                 if createPlots && nargin == 3
                     complete_state = state_LDF.filteredState;
                     
                     plotHandles(1) = figure;
                     p1 = subplot(2,1,1); hold on;
                     plot(complete_state(:,3),complete_state(:,6),'m.','MarkerSize',10);
                     ylabel('V_{gy} (m/s)', 'FontSize', 16);
                     set(gca, 'FontSize', 18); grid on;
                     
                     
                     p2 = subplot(2,1,2); hold on;
                     plot(complete_state(:,3),complete_state(:,6)./complete_state(:,3),'m.','MarkerSize',10);
                     ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
                     xlabel('y (m)', 'FontSize', 16);
                     set(gca, 'FontSize', 18); grid on;
                     
                     linkaxes([p1, p2],'x');
                         
                     state_subset = object.data4rrefEstimate(1).state4rrefEstimate;
                     plotHandles(2) = figure;
                     p1 = subplot(2,1,1); hold on;
                     plot(state_subset(:,3),state_subset(:,6),'r.','MarkerSize',15);
                     ylabel('V_{gy} (m/s)', 'FontSize', 16);
                     set(gca, 'FontSize', 18); grid on;
                     
                     p2 = subplot(2,1,2); hold on;
                     plot(state_subset(:,3),state_subset(:,6)./state_subset(:,3),'r.','MarkerSize',15);
                     ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
                     xlabel('y (m)', 'FontSize', 16);
                     set(gca, 'FontSize', 18); grid on;
%                      ylim([object.rref-1, obj.rref+1]);
                     
                     linkaxes([p1, p2],'x');
                 end
                 for ct=1:length(object.data4rrefEstimate) % For each identified segment
                     obj = object.data4rrefEstimate(ct);                     
                     obj.model = fitlm(obj.state4rrefEstimate(:,3),obj.state4rrefEstimate(:,6),'Intercept',false);
                     obj.rref = obj.model.Coefficients.Estimate;
                     obj.vmean = mean(obj.state4rrefEstimate(:,6));  
                     obj.ymean = mean(obj.state4rrefEstimate(:,3));
                     obj.meanVbyy = mean(obj.state4rrefEstimate(:,6)./obj.state4rrefEstimate(:,3));
                     obj.Rsquared = obj.model.Rsquared.Ordinary;
                     
                     obj.dof_analytical = log(abs(obj.state4rrefEstimate(end,3)/obj.state4rrefEstimate(1,3)))/obj.rref;
                     obj.dof_actual = diff(obj.state4rrefEstimate([1, end],1));
    %                  b = robustfit(obj.state4rrefEstimate(:,3), obj.state4rrefEstimate(:,6));
    %                  c = robustfit(obj.state4rrefEstimate(:,3), obj.state4rrefEstimate(:,6)./obj.state4rrefEstimate(:,3));
    %                  b1 = robustfit(obj.state4rrefEstimate(:,6), obj.state4rrefEstimate(:,9),[],[],'off');
    %                  obj.rref = b(2);
    %                  obj.rref = mean(obj.state4rrefEstimate(:,6)./obj.state4rrefEstimate(:,3));

                     if createPlots && nargin == 3
                         
                         state_subset = obj.state4rrefEstimate;

                         figure(plotHandles(1))
                         subplot(2,1,1);
                         plot(state_subset(:,3),state_subset(:,6),'r.','MarkerSize',10);
                         plot(state_subset(:,3),obj.rref*state_subset(:,3),'b','LineWidth',2);
                         
                         subplot(2,1,2);
                         plot(state_subset(:,3),state_subset(:,6)./state_subset(:,3),'r.','MarkerSize',10);
    %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
                         plot(state_subset([1,end],3),[obj.rref, obj.rref],'b','LineWidth',2);
                         

                         figure(plotHandles(2))
                         subplot(2,1,1);
                         plot(state_subset(:,3),state_subset(:,6),'r.','MarkerSize',15);
                         plot(state_subset(:,3),obj.rref*state_subset(:,3),'b','LineWidth',2);
                         
                         subplot(2,1,2);
                         plot(state_subset(:,3),state_subset(:,6)./state_subset(:,3),'r.','MarkerSize',15);
    %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
                         plot(state_subset([1,end],3),[obj.rref, obj.rref],'b','LineWidth',2);
                         


    %                      ylim([-7 1]);

    %                      p2 = subplot(2,1,2); hold on;
    %                      plot(complete_state(:,3),complete_state(:,9)./complete_state(:,6),'m','LineWidth',2);
    %                      ylabel('a_{y} / V_{gy} (1/s)', 'FontSize', 16);
    %                      xlabel('y (m)', 'FontSize', 16);
    %                      set(gca, 'FontSize', 18); grid on;
    %                      ylim([-40 40]);

                         
                     else
                         plotHandles = [];
                     end
                 end
                 
                 if createPlots && nargin == 3
                     figure(plotHandles(1))
                     subplot(2,1,1);
                     title(['R^2 = ' num2str([object.data4rrefEstimate.Rsquared],3)], 'FontSize', 14);
                     
                     figure(plotHandles(2))
                     subplot(2,1,1);
                     title(['R^2 = ' num2str([object.data4rrefEstimate.Rsquared],3)], 'FontSize', 14);
                 end
                 
             else
                 warning('First find data4rrefEstimate variable using GUI data');
                 plotHandles = [];
             end
             
         end
     end
end