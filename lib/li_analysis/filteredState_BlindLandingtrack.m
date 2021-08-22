classdef filteredState_BlindLandingtrack < handle
    properties (Access = public)
        filteredState = [];  % Butterworth filtered state - 10 columns (timestamp	x	y	z	xvel	yvel	zvel xacc yacc zacc)
        landingSide = '';
        
        % to store parameters for r* segments (these r* segments are
        % automatically computed
        rrefSegments = data4rrefEstimate.empty;
        
        rrefEntrySegments = data4rrefEntry.empty;
        
        decSegments = trackSegment.empty; % deceleration segments
        
%         rrefSegments_tw_fac = rrefEstimateInterval.empty; % For each track excerpt, storing rref segments for each factor and each time window
%         rrefSegments_best = rrefEstimateInterval.empty; % For each track excerpt, storing rref segments best among all factors and time windows
        
    end

    methods
        function obj = filteredState_BlindLandingtrack(filteredState, landingSide)
            assert(size(filteredState,2) == 10);
            obj.filteredState = filteredState;
            obj.landingSide = landingSide;
        end
        
        function plotHandles = plot_states(obj)
            % Plots V vs y and V/y vs y            
            
            
            complete_state = obj.filteredState(1:end,:);
            
            plotHandles = figure;
            p1 = subplot(2,1,1); hold on;
            plot(-complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('V (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            
            p2 = subplot(2,1,2); hold on;
            plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('r (1/s)', 'FontSize', 16);
            xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            
            linkaxes([p1, p2],'x');
            
        end
        
        function has = hasTakeoff(obj, landingDiscs)
            % Determines if the landing track has a non-takeoff start
            
            % landingDiscs - 1x2 LandingDisc array
            
            
            [~, ymin_indx] = min(obj.filteredState(:,3));
            state0 = obj.filteredState(ymin_indx,:);
            
            oppositeDisc = landingDiscs(arrayfun(@(x) ~strcmpi(x.side, obj.landingSide), landingDiscs));
            thisDisc = landingDiscs(arrayfun(@(x) strcmpi(x.side, obj.landingSide), landingDiscs));
            
            center = (oppositeDisc.center - thisDisc.center)';
            if strcmpi(thisDisc.side,'feeder')
                center(2) = -center(:, 2);
            end
            height_cylinder = 0.05;
            if any(filteredState_BlindLandingtrack.IsInsideCylinder(center - [0 height_cylinder/2 0], center + [0 height_cylinder/2 0], oppositeDisc.radius, obj.filteredState(:,2:4))) || ...
                    any(obj.filteredState(:,4) < -0.20)
                has = true;
            else
                has = false;
            end
        end
                
        function compute_rref(obj, params, factors, time_window)
            % params - [sigma_{rmean-c_{r vs y, 0-1}}, sigma_{m_{r vs y, 0-1}}] = [sigma1 sigma2]
            % factors - vector of some size
            % time_window - 2 by 1 vector = [min_gap max_gap]
            
            % For each point, intervals are found by looking ahead
            % min_gap:1:max_gap points
            
            % Within each such interval 3 lines are fit
            % First line between whole interval: r = m_{0-1} y + c_{0-1}
            % Second line between first half of the interval: r = m_{0-0.5} y + c_{0-0.5}
            % Third line between second half of the interval: r = m_{0.5-1} y + c_{0.5-1}
            % Time windows/intervals which satisfy the following constraint are
            % classsified as intervals with constant r (or simply r*).
            % Constraint: 
            
            assert(numel(params)==2 && numel(time_window)==2);
            rmse_func = @(idx, indices, data) sqrt(sum((data(indices(idx,1):indices(idx,2)) - mean(data(indices(idx,1):indices(idx,2)))).^2) ...
                                                   /(indices(idx,2)-indices(idx,1)+1));
            
            min_gap = time_window(1); % the first straight line is between current point and (current point + min_gap)th point
            max_gap = time_window(2);
            
%             % Taking track up to the first point before y = 2cm
%             indx = find(obj.filteredState(1:end,3)>-0.02, 1, 'first') - 1;
%             indx = min([indx size(obj.filteredState,1)]);
            indx = size(obj.filteredState,1);
            
            t_all = obj.filteredState(1:indx,1)-obj.filteredState(1,1);
            y_all = obj.filteredState(1:indx,3);
            v_all = obj.filteredState(1:indx,6);
            a_all = obj.filteredState(1:indx,9);
            r_all = v_all./y_all;
%             rdot_all = abs(a_all./y_all-(v_all./y_all).^2);

            x_all = obj.filteredState(1:indx,2);
            z_all = obj.filteredState(1:indx,4);
%             figure;
%             plot(t_all,r_all,'o');
%             text(t_all,r_all,string(1:length(r_all)));
%             title('t vs r')
            
%             figure;
%             plot(y_all,r_all,'o');
%             text(y_all,r_all,string(1:length(r_all)));
%             title('y vs r')
            
%             figure;
%             plot(t_all,rdot_all,'ro');
%             text(t_all,rdot_all,string(1:length(rdot_all)));
%             title('t vs rdot');
            
            
            for ct_factor=1:length(factors)
                obj.rrefSegments(ct_factor) = data4rrefEstimate();
                
                % Saving some parameters
                obj.rrefSegments(ct_factor).factor = factors(ct_factor);
                obj.rrefSegments(ct_factor).params = params;
            end
            
            N = indx;
            possible_indices = cell(length(factors),1); % possible intervals for each factor
            
            for ct=1:N-min_gap % for each point as a starting point for different lines through origin
                if v_all(ct)/y_all(ct) >= 0
                    continue;
                end
                
                end_point = min([ct+max_gap, N]);
                
                % Useful column vectors
%                 t = t_all(ct:end_point)-t_all(1);
                y = y_all(ct:end_point);
                v = v_all(ct:end_point);
                a = a_all(ct:end_point);
                r = v./y;
%                 rdot = abs(a./y-(v./y).^2);
                
                last_index = find(r > 0, 1) - 1; % last but one point where r becomes positive
                if isempty(last_index)
                    last_index = length(r);
                end
                
                coeff = zeros(last_index-min_gap,2); % [intercept slope]
                coeff_firsthalf = zeros(last_index-min_gap,2); 
                coeff_secondhalf = zeros(last_index-min_gap,2); 
                r_mean = zeros(last_index-min_gap,1); % Mean of r within each interval
                a_mean = zeros(last_index-min_gap,1); % Mean of a within each interval
                a_mean_firsthalf = zeros(last_index-min_gap,1); % mean of the ay in the first half of the interval
                a_mean_secondhalf = zeros(last_index-min_gap,1); % mean of the ay in the second half of the interval
%                 rdot_mean = zeros(last_index-min_gap,1); % Mean of r within each interval
%                 rdot_mean_firsthalf = zeros(last_index-min_gap,1); % Mean of r within each interval
%                 rdot_mean_secondhalf = zeros(last_index-min_gap,1); % Mean of r within each interval
                for ct1=min_gap+1:last_index
                    coeff(ct1-min_gap,:) = [ones(ct1,1) y(1:ct1)]\r(1:ct1);
                    r_mean(ct1-min_gap) = mean(r(1:ct1));
                    a_mean(ct1-min_gap) = mean(a(1:ct1));
%                     rdot_mean(ct1-min_gap) = mean(rdot(1:ct1));
                    if rem(ct1,2) == 0
%                         rdot_mean_firsthalf(ct1-min_gap) = mean(rdot(1:ct1/2));
%                         rdot_mean_secondhalf(ct1-min_gap) = mean(rdot(ct1/2+1:ct1));

                        a_mean_firsthalf(ct1-min_gap) = mean(a(1:ct1/2));
                        a_mean_secondhalf(ct1-min_gap) = mean(a(ct1/2+1:ct1));
                        
                        coeff_firsthalf(ct1-min_gap,:) = [ones(ct1/2,1) y(1:ct1/2)]\r(1:ct1/2);
                        coeff_secondhalf(ct1-min_gap,:) = [ones(ct1/2,1) y(ct1/2+1:ct1)]\r(ct1/2+1:ct1);
                    else
%                         rdot_mean_firsthalf(ct1-min_gap) = mean(rdot(1:floor(ct1/2)+1));
%                         rdot_mean_secondhalf(ct1-min_gap) = mean(rdot(floor(ct1/2)+1:ct1));
                        a_mean_firsthalf(ct1-min_gap) = mean(a(1:floor(ct1/2)+1));
                        a_mean_secondhalf(ct1-min_gap) = mean(a(floor(ct1/2)+1:ct1));  

                        coeff_firsthalf(ct1-min_gap,:) = [ones(floor(ct1/2)+1,1) y(1:floor(ct1/2)+1)]\r(1:floor(ct1/2)+1);
                        coeff_secondhalf(ct1-min_gap,:) = [ones(floor(ct1/2)+1,1) y(floor(ct1/2)+1:ct1)]\r(floor(ct1/2)+1:ct1); 
                    end
                end
                
                for ct_factor=1:length(factors)
                    factor = factors(ct_factor);
                    
%                     indx = find(abs(coeff(:,1)-r_mean) <= factor*params(1) & abs(coeff(:,2)) <= factor*params(2) & ...
%                                 abs(coeff_firsthalf(:,1)-r_mean) <= factor*2*params(1) & abs(coeff_firsthalf(:,2)) <= factor*2*params(2) & ...
%                                 abs(coeff_secondhalf(:,1)-r_mean) <= factor*2*params(1) & abs(coeff_secondhalf(:,2)) <= factor*2*params(2) & ...
%                                 a_mean <= 0 & a_mean >= -4 & ...
%                                 a_mean_firsthalf <= 0 & a_mean_firsthalf >= -4 & ...
%                                 a_mean_secondhalf <= 0 & a_mean_secondhalf >= -4);
                            
                    indx = find(abs(coeff(:,1)-r_mean) <= factor*params(1) & abs(coeff(:,2)) <= factor*params(2) & ...
                                abs(coeff_firsthalf(:,1)-r_mean) <= factor*2*params(1) & abs(coeff_firsthalf(:,2)) <= factor*2*params(2) & ...
                                abs(coeff_secondhalf(:,1)-r_mean) <= factor*2*params(1) & abs(coeff_secondhalf(:,2)) <= factor*2*params(2) & ...
                                a_mean <= 0 & ...
                                a_mean_firsthalf <= 0 & ...
                                a_mean_secondhalf <= 0);
                    if ~isempty(indx)
                        possible_indices{ct_factor} = [possible_indices{ct_factor}; ct*ones(numel(indx),1) ct+min_gap+indx-1];
                    end
                    
                end
                
            end
            
            for ct_factor=1:length(factors)
                rmse_constFit = arrayfun(@(idx) rmse_func(idx, possible_indices{ct_factor}, r_all), 1:size(possible_indices{ct_factor},1));
                if ~isempty(possible_indices{ct_factor})
                    intervalArray = Interval.createIntervalArray(possible_indices{ct_factor}, rmse_constFit');
                    
                    
                    
                    % Time-window independent best intervals
                    chosenIntervals_ti = Interval.findRepresentativeIntervals(intervalArray);
                    intervals_ti = possible_indices{ct_factor}(chosenIntervals_ti,:);
                    obj.rrefSegments(ct_factor).intervals_ti = intervals_ti;
                    obj.rrefSegments(ct_factor).rref_ti = arrayfun(@(i,j) sum(y_all(i:j).*v_all(i:j))/sum(y_all(i:j).^2),intervals_ti(:,1),intervals_ti(:,2));
                    obj.rrefSegments(ct_factor).ymean_ti = arrayfun(@(i,j) mean(y_all(i:j)),intervals_ti(:,1),intervals_ti(:,2));
                    obj.rrefSegments(ct_factor).vmean_ti = arrayfun(@(i,j) mean(v_all(i:j)),intervals_ti(:,1),intervals_ti(:,2));
                    obj.rrefSegments(ct_factor).rmean_ti = arrayfun(@(i,j) mean(r_all(i:j)),intervals_ti(:,1),intervals_ti(:,2));
                    coeff = arrayfun(@(i,j) [[ones(j-i+1,1) y_all(i:j)]\r_all(i:j)]',intervals_ti(:,1),intervals_ti(:,2),'UniformOutput',false);
                    coeff = vertcat(coeff{:});
                    obj.rrefSegments(ct_factor).const_rvsy_ti = coeff(:,1);
                    obj.rrefSegments(ct_factor).slope_rvsy_ti = coeff(:,2);
                    
                    obj.rrefSegments(ct_factor).xTravelled_ti = arrayfun(@(i,j) abs(max(x_all(i:j)) - min(x_all(i:j))), intervals_ti(:,1), intervals_ti(:,2)) ;
                    obj.rrefSegments(ct_factor).yTravelled_ti = arrayfun(@(i,j) abs(max(y_all(i:j)) - min(y_all(i:j))), intervals_ti(:,1), intervals_ti(:,2)) ;
                    obj.rrefSegments(ct_factor).zTravelled_ti = arrayfun(@(i,j) abs(max(z_all(i:j)) - min(z_all(i:j))), intervals_ti(:,1), intervals_ti(:,2)) ;
                    obj.rrefSegments(ct_factor).xmean_ti = arrayfun(@(i,j) mean(x_all(i:j)), intervals_ti(:,1), intervals_ti(:,2)) ;
                    obj.rrefSegments(ct_factor).zmean_ti = arrayfun(@(i,j) mean(z_all(i:j)), intervals_ti(:,1), intervals_ti(:,2)) ;
                    obj.rrefSegments(ct_factor).yrange = arrayfun(@(i,j) abs(y_all(i) - y_all(j)), intervals_ti(:,1), intervals_ti(:,2)) ;
                    obj.rrefSegments(ct_factor).ymid = arrayfun(@(i,j) y_all(i) + abs(y_all(i) - y_all(j))/2, intervals_ti(:,1), intervals_ti(:,2)) ;
                    obj.rrefSegments(ct_factor).fd_analytical_ti = arrayfun(@(i,j) log(abs(y_all(j)/y_all(i))), intervals_ti(:,1), intervals_ti(:,2))./obj.rrefSegments(ct_factor).rref_ti; % analytically computed flight duration 
                    obj.rrefSegments(ct_factor).fd_actual_ti = arrayfun(@(i,j) diff(t_all([i j])), intervals_ti(:,1), intervals_ti(:,2)); % actual flight duration
                    
                    % Time-window dependent best intervals for a particular
                    % factor
                    chosenIntervals_td = Interval.findRepresentativeIntervalsPerTimeWindow(intervalArray);
                    intervals_td = possible_indices{ct_factor}(chosenIntervals_td,:);
                    obj.rrefSegments(ct_factor).intervals_td = intervals_td;
                    obj.rrefSegments(ct_factor).rref_td = arrayfun(@(i,j) sum(y_all(i:j).*v_all(i:j))/sum(y_all(i:j).^2),intervals_td(:,1),intervals_td(:,2));
                    obj.rrefSegments(ct_factor).ymean_td = arrayfun(@(i,j) mean(y_all(i:j)),intervals_td(:,1),intervals_td(:,2));
                    obj.rrefSegments(ct_factor).vmean_td = arrayfun(@(i,j) mean(v_all(i:j)),intervals_td(:,1),intervals_td(:,2));
%                     
                    
                end
            end
            
            % Throw all intervals (different factors and time-windows) together and find best representatives
            
            
        end
        
        function plotHandles = plot_rrefs(obj, factor)
            % Plots V vs y and V/y vs y highling change in r* within the same track
            % excerpt with positive V and positive y
            % This function plots time-window independent best intervals
            
            % factor - for which threshold factor plots are required to be
            % produced
            
            assert(~isempty(obj.rrefSegments) && length(factor) == 1);
            ct_factor = find(abs([obj.rrefSegments.factor] - factor) < 1e-6);
            if isempty(ct_factor)
                error('Can NOT find the r* intervals for the asked factor.');
            end
            intervals = obj.rrefSegments(ct_factor).intervals_ti;
            if isempty(intervals)
                plotHandles = [];
                return;
            end
            
            dataPerTrackExcerpt = struct.empty;
            for ct=1:size(intervals,1)
                dataPerTrackExcerpt(ct).state4rrefEstimate = obj.filteredState(intervals(ct,1):intervals(ct,2),:);
                dataPerTrackExcerpt(ct).rref = obj.rrefSegments(ct_factor).rref_ti(ct);%mean(dataPerTrackExcerpt(ct).state4rrefEstimate(:,6)./dataPerTrackExcerpt(ct).state4rrefEstimate(:,3));
            end
            data4rref = vertcat(dataPerTrackExcerpt.state4rrefEstimate);
            
            ymin = min(data4rref(:,3));
            [~, indx_ymin] = min(abs(obj.filteredState(:,3)-ymin));
            time_ymin = obj.filteredState(indx_ymin,1);
            % Start from 0.06s data earlier than dataTrackExcerpt.filteredState(indx_ymax,1)
            indx_start = find(obj.filteredState(:,1)-time_ymin<-0.06,1);
            if isempty(indx_start)
                indx_start = 1;
            end
            complete_state = obj.filteredState(indx_start:end,:);
            
            plotHandles = figure;
            p1 = subplot(2,1,1); hold on;
            plot(-complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('V (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3)], 'FontSize', 16);
            %     ylim([-0.2 p1.YLim(2)]);
            %     ylim([-0.2 0.6]);
            %     yticks([-0.2:0.2:0.6]);
            %     xticks([p1.XLim(1):0.1:0]);
            
            
            p2 = subplot(2,1,2); hold on;
            %     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('r (1/s)', 'FontSize', 16);
            %     ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
            xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            ylim([0 ceil(-1*ceil(min(data4rref(:,6)./data4rref(:,3))-1))]);
            yticks(0:1:ceil(-1*ceil(min(data4rref(:,6)./data4rref(:,3))-1)));
            %     xticks(p1.XLim(1):0.1:0);
            
            linkaxes([p1, p2],'x');
            if abs(floor((min(data4rref(:,3))/0.05))*0.05 - min(data4rref(:,3))) > 0.025
                xmin = floor((min(data4rref(:,3))/0.05))*0.05;
            else
                xmin = floor((min(data4rref(:,3))/0.05))*0.05-0.05;
            end
            xlim([0 -xmin]);
            %     xlim([0 Inf]);
            
            for ct=1:size(intervals,1)
                obj = dataPerTrackExcerpt(ct);
                state_subset = obj.state4rrefEstimate;
                figure(plotHandles(1))
                subplot(2,1,1);
                plot(-state_subset(:,3),state_subset(:,6),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(-state_subset(:,3),-obj.rref*-state_subset(:,3),'LineWidth',2,'Color',[69 117 180]./255);
                plot([0 -state_subset(end,3)],[0 -obj.rref*-state_subset(end,3)],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                subplot(2,1,2);
                plot(-state_subset(:,3),state_subset(:,6)./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
                plot(-state_subset([1,end],3),[-obj.rref, -obj.rref],'LineWidth',2,'Color',[69 117 180]./255);
                plot([p1.XLim(1) -state_subset(1,3)],[-obj.rref -obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
            end
            
        end
        
        function plotHandles = plot_rrefs_with_acc(obj, factor)
            % Plots V vs y, V/y vs y and a vs y highling change in r* within the same track
            % excerpt with positive V and positive y
            % This function plots time-window independent best intervals
            
            % factor - for which threshold factor plots are required to be
            % produced
            
            assert(~isempty(obj.rrefSegments) && length(factor) == 1);
            ct_factor = find(abs([obj.rrefSegments.factor] - factor) < 1e-6);
            if isempty(ct_factor)
                error('Can NOT find the r* intervals for the asked factor.');
            end
            intervals = obj.rrefSegments(ct_factor).intervals_ti;
            if isempty(intervals)
                plotHandles = [];
                return;
            end
            
            dataPerTrackExcerpt = struct.empty;
            for ct=1:size(intervals,1)
                dataPerTrackExcerpt(ct).state4rrefEstimate = obj.filteredState(intervals(ct,1):intervals(ct,2),:);
                dataPerTrackExcerpt(ct).rref = obj.rrefSegments(ct_factor).rref_ti(ct);%mean(dataPerTrackExcerpt(ct).state4rrefEstimate(:,6)./dataPerTrackExcerpt(ct).state4rrefEstimate(:,3));
            end
            data4rref = vertcat(dataPerTrackExcerpt.state4rrefEstimate);
            
            ymin = min(data4rref(:,3));
            [~, indx_ymin] = min(abs(obj.filteredState(:,3)-ymin));
            time_ymin = obj.filteredState(indx_ymin,1);
            % Start from 0.06s data earlier than dataTrackExcerpt.filteredState(indx_ymax,1)
            indx_start = find(obj.filteredState(:,1)-time_ymin<-0.06,1);
            if isempty(indx_start)
                indx_start = 1;
            end
            complete_state = obj.filteredState(indx_start:end,:);
            
            plotHandles = figure;
            p1 = subplot(3,1,1); hold on;
            plot(-complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('V (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3)], 'FontSize', 16);
            %     ylim([-0.2 p1.YLim(2)]);
            %     ylim([-0.2 0.6]);
            %     yticks([-0.2:0.2:0.6]);
            %     xticks([p1.XLim(1):0.1:0]);
            
            
            p2 = subplot(3,1,2); hold on;
            %     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('r (1/s)', 'FontSize', 16);
            %     ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
            xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            ylim([0 ceil(-1*ceil(min(data4rref(:,6)./data4rref(:,3))-1))]);
            yticks(0:1:ceil(-1*ceil(min(data4rref(:,6)./data4rref(:,3))-1)));
            %     xticks(p1.XLim(1):0.1:0);
            
            p3 = subplot(3,1,3); hold on;
            %     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-complete_state(:,3),complete_state(:,9),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('a (ms^{-2})', 'FontSize', 16);
            xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
%             ylim([0 ceil(-1*ceil(min(data4rref(:,6)./data4rref(:,3))-1))]);
%             yticks(0:1:ceil(-1*ceil(min(data4rref(:,6)./data4rref(:,3))-1)));
            
            linkaxes([p1, p2, p3],'x');
            if abs(floor((min(data4rref(:,3))/0.05))*0.05 - min(data4rref(:,3))) > 0.025
                xmin = floor((min(data4rref(:,3))/0.05))*0.05;
            else
                xmin = floor((min(data4rref(:,3))/0.05))*0.05-0.05;
            end
            xlim([0 -xmin]);
            %     xlim([0 Inf]);
            
            for ct=1:size(intervals,1)
                obj = dataPerTrackExcerpt(ct);
                state_subset = obj.state4rrefEstimate;
                figure(plotHandles(1))
                subplot(3,1,1);
                plot(-state_subset(:,3),state_subset(:,6),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(-state_subset(:,3),-obj.rref*-state_subset(:,3),'LineWidth',2,'Color',[69 117 180]./255);
                plot([0 -state_subset(end,3)],[0 -obj.rref*-state_subset(end,3)],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                subplot(3,1,2);
                plot(-state_subset(:,3),state_subset(:,6)./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
                plot(-state_subset([1,end],3),[-obj.rref, -obj.rref],'LineWidth',2,'Color',[69 117 180]./255);
                plot([p1.XLim(1) -state_subset(1,3)],[-obj.rref -obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
            end
            
        end
        
        function compute_rref_with3dspeed(obj, params, factors, time_window)
            % params - [sigma_{rmean-c_{r vs y, 0-1}}, sigma_{m_{r vs y, 0-1}}] = [sigma1 sigma2]
            % factors - vector of some size
            % time_window - 2 by 1 vector = [min_gap max_gap]
            
            % For each point, intervals are found by looking ahead
            % min_gap:1:max_gap points
            
            % Within each such interval 3 lines are fit
            % First line between whole interval: r = m_{0-1} y + c_{0-1}
            % Second line between first half of the interval: r = m_{0-0.5} y + c_{0-0.5}
            % Third line between second half of the interval: r = m_{0.5-1} y + c_{0.5-1}
            % Time windows/intervals which satisfy the following constraint are
            % classsified as intervals with constant r (or simply r*).
            % Constraint: 
            
            assert(numel(params)==2 && numel(time_window)==2);
            rmse_func = @(idx, indices, data) sqrt(sum((data(indices(idx,1):indices(idx,2)) - mean(data(indices(idx,1):indices(idx,2)))).^2) ...
                                                   /(indices(idx,2)-indices(idx,1)+1));
            
            min_gap = time_window(1); % the first straight line is between current point and (current point + min_gap)th point
            max_gap = time_window(2);
            
%             % Taking track up to the first point before y = 2cm
%             indx = find(obj.filteredState(1:end,3)>-0.02, 1, 'first') - 1;
%             indx = min([indx size(obj.filteredState,1)]);
            indx = size(obj.filteredState,1);
            
            t_all = obj.filteredState(1:indx,1)-obj.filteredState(1,1);
            y_all = obj.filteredState(1:indx,3);
            vy_all = obj.filteredState(1:indx,6);
            v_all = (sum(obj.filteredState(1:indx,5:7).^2,2)).^0.5;
            a_all = obj.filteredState(1:indx,9);
            r_all = v_all./y_all;
%             rdot_all = abs(a_all./y_all-(v_all./y_all).^2);

            x_all = obj.filteredState(1:indx,2);
            z_all = obj.filteredState(1:indx,4);
%             figure;
%             plot(t_all,r_all,'o');
%             text(t_all,r_all,string(1:length(r_all)));
%             title('t vs r')
            
%             figure;
%             plot(y_all,r_all,'o');
%             text(y_all,r_all,string(1:length(r_all)));
%             title('y vs r')
            
%             figure;
%             plot(t_all,rdot_all,'ro');
%             text(t_all,rdot_all,string(1:length(rdot_all)));
%             title('t vs rdot');
            
            
            for ct_factor=1:length(factors)
                obj.rrefSegments(ct_factor) = data4rrefEstimate();
                
                % Saving some parameters
                obj.rrefSegments(ct_factor).factor = factors(ct_factor);
                obj.rrefSegments(ct_factor).params = params;
            end
            
            N = indx;
            possible_indices = cell(length(factors),1); % possible intervals for each factor
            
            for ct=1:N-min_gap % for each point as a starting point for different lines through origin
                if vy_all(ct)/y_all(ct) >= 0
                    continue;
                end
                
                end_point = min([ct+max_gap, N]);
                
                % Useful column vectors
%                 t = t_all(ct:end_point)-t_all(1);
                y = y_all(ct:end_point);
                v = v_all(ct:end_point);
                a = a_all(ct:end_point);
                r = v./y;
%                 rdot = abs(a./y-(v./y).^2);
                
                last_index = find(vy_all(ct:end_point)./y_all(ct:end_point) > 0, 1) - 1; % last but one point where r becomes positive
                if isempty(last_index)
                    last_index = length(r);
                end
                
                coeff = zeros(last_index-min_gap,2); % [intercept slope]
                coeff_firsthalf = zeros(last_index-min_gap,2); 
                coeff_secondhalf = zeros(last_index-min_gap,2); 
                r_mean = zeros(last_index-min_gap,1); % Mean of r within each interval
                a_mean = zeros(last_index-min_gap,1); % Mean of a within each interval
                a_mean_firsthalf = zeros(last_index-min_gap,1); % mean of the ay in the first half of the interval
                a_mean_secondhalf = zeros(last_index-min_gap,1); % mean of the ay in the second half of the interval
%                 rdot_mean = zeros(last_index-min_gap,1); % Mean of r within each interval
%                 rdot_mean_firsthalf = zeros(last_index-min_gap,1); % Mean of r within each interval
%                 rdot_mean_secondhalf = zeros(last_index-min_gap,1); % Mean of r within each interval
                for ct1=min_gap+1:last_index
                    coeff(ct1-min_gap,:) = [ones(ct1,1) y(1:ct1)]\r(1:ct1);
                    r_mean(ct1-min_gap) = mean(r(1:ct1));
                    a_mean(ct1-min_gap) = mean(a(1:ct1));
%                     rdot_mean(ct1-min_gap) = mean(rdot(1:ct1));
                    if rem(ct1,2) == 0
%                         rdot_mean_firsthalf(ct1-min_gap) = mean(rdot(1:ct1/2));
%                         rdot_mean_secondhalf(ct1-min_gap) = mean(rdot(ct1/2+1:ct1));

                        a_mean_firsthalf(ct1-min_gap) = mean(a(1:ct1/2));
                        a_mean_secondhalf(ct1-min_gap) = mean(a(ct1/2+1:ct1));
                        
                        coeff_firsthalf(ct1-min_gap,:) = [ones(ct1/2,1) y(1:ct1/2)]\r(1:ct1/2);
                        coeff_secondhalf(ct1-min_gap,:) = [ones(ct1/2,1) y(ct1/2+1:ct1)]\r(ct1/2+1:ct1);
                    else
%                         rdot_mean_firsthalf(ct1-min_gap) = mean(rdot(1:floor(ct1/2)+1));
%                         rdot_mean_secondhalf(ct1-min_gap) = mean(rdot(floor(ct1/2)+1:ct1));
                        a_mean_firsthalf(ct1-min_gap) = mean(a(1:floor(ct1/2)+1));
                        a_mean_secondhalf(ct1-min_gap) = mean(a(floor(ct1/2)+1:ct1));  

                        coeff_firsthalf(ct1-min_gap,:) = [ones(floor(ct1/2)+1,1) y(1:floor(ct1/2)+1)]\r(1:floor(ct1/2)+1);
                        coeff_secondhalf(ct1-min_gap,:) = [ones(floor(ct1/2)+1,1) y(floor(ct1/2)+1:ct1)]\r(floor(ct1/2)+1:ct1); 
                    end
                end
                
                for ct_factor=1:length(factors)
                    factor = factors(ct_factor);
                    
%                     indx = find(abs(coeff(:,1)-r_mean) <= factor*params(1) & abs(coeff(:,2)) <= factor*params(2) & ...
%                                 abs(coeff_firsthalf(:,1)-r_mean) <= factor*2*params(1) & abs(coeff_firsthalf(:,2)) <= factor*2*params(2) & ...
%                                 abs(coeff_secondhalf(:,1)-r_mean) <= factor*2*params(1) & abs(coeff_secondhalf(:,2)) <= factor*2*params(2) & ...
%                                 a_mean <= 0 & a_mean >= -4 & ...
%                                 a_mean_firsthalf <= 0 & a_mean_firsthalf >= -4 & ...
%                                 a_mean_secondhalf <= 0 & a_mean_secondhalf >= -4);
                            
                    indx = find(abs(coeff(:,1)-r_mean) <= factor*params(1) & abs(coeff(:,2)) <= factor*params(2) & ...
                                abs(coeff_firsthalf(:,1)-r_mean) <= factor*2*params(1) & abs(coeff_firsthalf(:,2)) <= factor*2*params(2) & ...
                                abs(coeff_secondhalf(:,1)-r_mean) <= factor*2*params(1) & abs(coeff_secondhalf(:,2)) <= factor*2*params(2) & ...
                                a_mean <= 0 & ...
                                a_mean_firsthalf <= 0 & ...
                                a_mean_secondhalf <= 0);
                    if ~isempty(indx)
                        possible_indices{ct_factor} = [possible_indices{ct_factor}; ct*ones(numel(indx),1) ct+min_gap+indx-1];
                    end
                    
                end
                
            end
            
            for ct_factor=1:length(factors)
                rmse_constFit = arrayfun(@(idx) rmse_func(idx, possible_indices{ct_factor}, r_all), 1:size(possible_indices{ct_factor},1));
                if ~isempty(possible_indices{ct_factor})
                    intervalArray = Interval.createIntervalArray(possible_indices{ct_factor}, rmse_constFit');
                    
                    
                    
                    % Time-window independent best intervals
                    chosenIntervals_ti = Interval.findRepresentativeIntervals(intervalArray);
                    intervals_ti = possible_indices{ct_factor}(chosenIntervals_ti,:);
                    obj.rrefSegments(ct_factor).intervals_ti = intervals_ti;
                    obj.rrefSegments(ct_factor).rref_ti = arrayfun(@(i,j) sum(y_all(i:j).*v_all(i:j))/sum(y_all(i:j).^2),intervals_ti(:,1),intervals_ti(:,2));
                    obj.rrefSegments(ct_factor).ymean_ti = arrayfun(@(i,j) mean(y_all(i:j)),intervals_ti(:,1),intervals_ti(:,2));
                    obj.rrefSegments(ct_factor).vmean_ti = arrayfun(@(i,j) mean(v_all(i:j)),intervals_ti(:,1),intervals_ti(:,2));
                    obj.rrefSegments(ct_factor).rmean_ti = arrayfun(@(i,j) mean(r_all(i:j)),intervals_ti(:,1),intervals_ti(:,2));
                    coeff = arrayfun(@(i,j) [[ones(j-i+1,1) y_all(i:j)]\r_all(i:j)]',intervals_ti(:,1),intervals_ti(:,2),'UniformOutput',false);
                    coeff = vertcat(coeff{:});
                    obj.rrefSegments(ct_factor).const_rvsy_ti = coeff(:,1);
                    obj.rrefSegments(ct_factor).slope_rvsy_ti = coeff(:,2);
                    
                    obj.rrefSegments(ct_factor).xTravelled_ti = arrayfun(@(i,j) abs(max(x_all(i:j)) - min(x_all(i:j))), intervals_ti(:,1), intervals_ti(:,2)) ;
                    obj.rrefSegments(ct_factor).yTravelled_ti = arrayfun(@(i,j) abs(max(y_all(i:j)) - min(y_all(i:j))), intervals_ti(:,1), intervals_ti(:,2)) ;
                    obj.rrefSegments(ct_factor).zTravelled_ti = arrayfun(@(i,j) abs(max(z_all(i:j)) - min(z_all(i:j))), intervals_ti(:,1), intervals_ti(:,2)) ;
                    obj.rrefSegments(ct_factor).xmean_ti = arrayfun(@(i,j) mean(x_all(i:j)), intervals_ti(:,1), intervals_ti(:,2)) ;
                    obj.rrefSegments(ct_factor).zmean_ti = arrayfun(@(i,j) mean(z_all(i:j)), intervals_ti(:,1), intervals_ti(:,2)) ;
                    obj.rrefSegments(ct_factor).fd_analytical_ti = arrayfun(@(i,j) log(abs(y_all(j)/y_all(i))), intervals_ti(:,1), intervals_ti(:,2))./obj.rrefSegments(ct_factor).rref_ti; % analytically computed flight duration 
                    obj.rrefSegments(ct_factor).fd_actual_ti = arrayfun(@(i,j) diff(t_all([i j])), intervals_ti(:,1), intervals_ti(:,2)); % actual flight duration
                    
                    % Time-window dependent best intervals for a particular
                    % factor
                    chosenIntervals_td = Interval.findRepresentativeIntervalsPerTimeWindow(intervalArray);
                    intervals_td = possible_indices{ct_factor}(chosenIntervals_td,:);
                    obj.rrefSegments(ct_factor).intervals_td = intervals_td;
                    obj.rrefSegments(ct_factor).rref_td = arrayfun(@(i,j) sum(y_all(i:j).*v_all(i:j))/sum(y_all(i:j).^2),intervals_td(:,1),intervals_td(:,2));
                    obj.rrefSegments(ct_factor).ymean_td = arrayfun(@(i,j) mean(y_all(i:j)),intervals_td(:,1),intervals_td(:,2));
                    obj.rrefSegments(ct_factor).vmean_td = arrayfun(@(i,j) mean(v_all(i:j)),intervals_td(:,1),intervals_td(:,2));
%                     
                    
                end
            end
            
            % Throw all intervals (different factors and time-windows) together and find best representatives
            
            
        end
        
        
        function plotHandles = plot_rrefs_with3dspeed(obj, factor)
            % Plots V vs y and V/y vs y highling change in r* within the same track
            % excerpt with positive V and positive y. Also plots 3d
            % velocity
            % This function plots time-window independent best intervals
            
            % factor - for which threshold factor plots are required to be
            % produced
            
            assert(~isempty(obj.rrefSegments) && length(factor) == 1);
            ct_factor = find(abs([obj.rrefSegments.factor] - factor) < 1e-6);
            if isempty(ct_factor)
                error('Can NOT find the r* intervals for the asked factor.');
            end
            intervals = obj.rrefSegments(ct_factor).intervals_ti;
            if isempty(intervals)
                plotHandles = [];
                return;
            end
            
            dataPerTrackExcerpt = struct.empty;
            for ct=1:size(intervals,1)
                dataPerTrackExcerpt(ct).state4rrefEstimate = obj.filteredState(intervals(ct,1):intervals(ct,2),:);
                dataPerTrackExcerpt(ct).rref = obj.rrefSegments(ct_factor).rref_ti(ct);%mean(dataPerTrackExcerpt(ct).state4rrefEstimate(:,6)./dataPerTrackExcerpt(ct).state4rrefEstimate(:,3));
            end
            data4rref = vertcat(dataPerTrackExcerpt.state4rrefEstimate);
            data4rref_speed3D = (sum(data4rref(:,[5 6 7]).^2,2)).^0.5;
            
            ymin = min(data4rref(:,3));
            [~, indx_ymin] = min(abs(obj.filteredState(:,3)-ymin));
            time_ymin = obj.filteredState(indx_ymin,1);
            % Start from 0.06s data earlier than dataTrackExcerpt.filteredState(indx_ymax,1)
            indx_start = find(obj.filteredState(:,1)-time_ymin<-0.06,1);
            if isempty(indx_start)
                indx_start = 1;
            end
            complete_state = obj.filteredState(indx_start:end,:);
            speed3D = (sum(complete_state(:,[5 6 7]).^2,2)).^0.5;
            
            plotHandles = figure;
            p1 = subplot(2,1,1); hold on;
            plot(-complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-complete_state(:,3),speed3D,'.','MarkerSize',10,'MarkerFaceColor',[170 255 195]./255, 'MarkerEdgeColor',[170 255 195]./255');
            ylabel('V (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3) ', 3D speed: green, Vy: orange'], 'FontSize', 16);
            %     ylim([-0.2 p1.YLim(2)]);
            %     ylim([-0.2 0.6]);
            %     yticks([-0.2:0.2:0.6]);
            %     xticks([p1.XLim(1):0.1:0]);
            
            
            p2 = subplot(2,1,2); hold on;
            %     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-complete_state(:,3),speed3D./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[170 255 195]./255, 'MarkerEdgeColor',[170 255 195]./255');
            ylabel('r (1/s)', 'FontSize', 16);
            %     ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
            xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            ylim([0 ceil(-1*ceil(min([data4rref(:,6)./data4rref(:,3); data4rref_speed3D./data4rref(:,3)])-1))]);
            yticks(0:1:ceil(-1*ceil(min(data4rref(:,6)./data4rref(:,3))-1)));
            %     xticks(p1.XLim(1):0.1:0);
            
            linkaxes([p1, p2],'x');
            if abs(floor((min(data4rref(:,3))/0.05))*0.05 - min(data4rref(:,3))) > 0.025
                xmin = floor((min(data4rref(:,3))/0.05))*0.05;
            else
                xmin = floor((min(data4rref(:,3))/0.05))*0.05-0.05;
            end
            xlim([0 -xmin]);
            %     xlim([0 Inf]);
            
            for ct=1:size(intervals,1)
                obj = dataPerTrackExcerpt(ct);
                state_subset = obj.state4rrefEstimate;
                speed3D = (sum(state_subset(:,[5 6 7]).^2,2)).^0.5;
                figure(plotHandles(1))
                subplot(2,1,1);
                plot(-state_subset(:,3),speed3D,'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(-state_subset(:,3),-obj.rref*-state_subset(:,3),'LineWidth',2,'Color',[69 117 180]./255);
                plot([0 -state_subset(end,3)],[0 -obj.rref*-state_subset(end,3)],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                subplot(2,1,2);
                plot(-state_subset(:,3),speed3D./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
                plot(-state_subset([1,end],3),[-obj.rref, -obj.rref],'LineWidth',2,'Color',[69 117 180]./255);
                plot([p1.XLim(1) -state_subset(1,3)],[-obj.rref -obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
            end
            
        end
        
        
        function setLandingSide(obj)
            % To set landing side for all extracted rref segments
            
            if strcmpi(obj.landingSide, 'Hive')
                [obj.rrefSegments.side] = deal(1);
                [obj.rrefEntrySegments.side] = deal(1);
            elseif strcmpi(obj.landingSide, 'Feeder')
                [obj.rrefSegments.side] = deal(2);
                [obj.rrefEntrySegments.side] = deal(2);
            end
        end
        
        function compute_dec_segs(obj)
            % Compute deceleration segments in a track
            
            % Within each such interval 3 lines are fit
            % First line between whole interval: r = m_{0-1} y + c_{0-1}
            % Second line between first half of the interval: r = m_{0-0.5} y + c_{0-0.5}
            % Third line between second half of the interval: r = m_{0.5-1} y + c_{0.5-1}
            % Time windows/intervals which satisfy the following constraint are
            % classsified as intervals with constant r (or simply r*).
            % Constraint: 
                        
            min_gap = 15; % the first straight line is between current point and (current point + min_gap)th point
            max_gap = 100; %size(obj.filteredState,1);
            
%             % Taking track up to the first point before y = 2cm
%             indx = find(obj.filteredState(1:end,3)>-0.02, 1, 'first') - 1;
%             indx = min([indx size(obj.filteredState,1)]);
            indx = size(obj.filteredState,1);
            
            t_all = obj.filteredState(1:indx,1)-obj.filteredState(1,1);
            y_all = obj.filteredState(1:indx,3);
            v_all = obj.filteredState(1:indx,6);
            a_all = obj.filteredState(1:indx,9);
            r_all = v_all./y_all;
%             rdot_all = abs(a_all./y_all-(v_all./y_all).^2);

            x_all = obj.filteredState(1:indx,2);
            z_all = obj.filteredState(1:indx,4);
%             figure;
%             plot(t_all,r_all,'o');
%             text(t_all,r_all,string(1:length(r_all)));
%             title('t vs r')
            
%             figure;
%             plot(y_all,r_all,'o');
%             text(y_all,r_all,string(1:length(r_all)));
%             title('y vs r')
            
%             figure;
%             plot(t_all,rdot_all,'ro');
%             text(t_all,rdot_all,string(1:length(rdot_all)));
%             title('t vs rdot');

            figure;
            plot(t_all,a_all,'ro');
            text(t_all,a_all,string(1:length(a_all)));
            title('t vs a');
            
            figure;
            plot(y_all,v_all,'ro');
            text(y_all,v_all,string(1:length(v_all)));
            title('y vs v');
            
            obj.decSegments = trackSegment();
            
            N = indx;
            possible_indices = []; % possible intervals 
            a_score = [];
            for ct=1:N-min_gap % for each point as a starting point for different lines through origin
                if v_all(ct)/y_all(ct) >= 0
                    continue;
                end
                
                end_point = min([ct+max_gap, N]);
                
                % Useful column vectors
%                 t = t_all(ct:end_point)-t_all(1);
                y = y_all(ct:end_point);
                v = v_all(ct:end_point);
                a = a_all(ct:end_point);
                r = v./y;
%                 rdot = abs(a./y-(v./y).^2);
                
                last_index = find(r > 0, 1) - 1; % last but one point where r becomes positive
                if isempty(last_index)
                    last_index = length(r);
                end
                
                isVpos = zeros(last_index-min_gap,1); % True if all V are positive
                
                a_mean = zeros(last_index-min_gap,1); % Mean of a within each interval
                a_mean_firsthalf = zeros(last_index-min_gap,1); % mean of the ay in the first half of the interval
                a_mean_secondhalf = zeros(last_index-min_gap,1); % mean of the ay in the second half of the interval

                for ct1=min_gap+1:last_index
                    isVpos(ct1-min_gap) = all(v(1:ct1)>0);
                    
                    a_mean(ct1-min_gap) = mean(a(1:ct1));
                    if rem(ct1,2) == 0
                        a_mean_firsthalf(ct1-min_gap) = mean(a(1:ct1/2));
                        a_mean_secondhalf(ct1-min_gap) = mean(a(ct1/2+1:ct1));
                    else
                        a_mean_firsthalf(ct1-min_gap) = mean(a(1:floor(ct1/2)+1));
                        a_mean_secondhalf(ct1-min_gap) = mean(a(floor(ct1/2)+1:ct1));  
                    end
                end
                
                indx = find(isVpos & ...
                    a_mean <= 0 & ...
                    a_mean_firsthalf <= 0 & ...
                    a_mean_secondhalf <= 0);
                
                if ~isempty(indx)
                    possible_indices = [possible_indices; ct*ones(numel(indx),1) ct+min_gap+indx-1];
                    dummy = [a_mean(indx) a_mean_firsthalf(indx) a_mean_secondhalf(indx)];
                    a_score = [a_score; min(dummy,[],2)];
                end
            end
            
            
            if ~isempty(possible_indices)
                len_inv = 1./(diff(possible_indices,[],2)+1);
                intervalArray = Interval.createIntervalArray(possible_indices, len_inv);
                
                % Time-window independent best intervals
                chosenIntervals_ti = Interval.findRepresentativeIntervals(intervalArray);
                intervals_ti = possible_indices(chosenIntervals_ti,:);
                obj.decSegments.intervals = intervals_ti;
                
            end
        end
        
        function compute_params_basedon_3dspeed(obj)
            for ct=1:length(obj.rrefSegments)
                intervals = obj.rrefSegments(ct).intervals_ti;
                obj.rrefSegments(ct).speed3d_mean_ti = arrayfun(@(x) mean((sum(obj.filteredState(intervals(x,1):intervals(x,2),5:7).^2, 2)).^0.5),1:size(intervals,1))';
                obj.rrefSegments(ct).rmean_speed3d_ti = arrayfun(@(x) mean((sum(obj.filteredState(intervals(x,1):intervals(x,2),5:7).^2, 2)).^0.5./obj.filteredState(intervals(x,1):intervals(x,2),3)),1:size(intervals,1))';
            end
        end
        
        function find_rrefEntry(obj)
            % Goes through obj.rrefSegments and finds track segments
            % corresponding to entry into constant r.
            % monotonic increase or decrease (decided by first point just
            % before start of constant-r) decides the start of entry
            
            if isempty(obj.rrefSegments)
                return
            end
            
            indx = size(obj.filteredState,1);
            t_all = obj.filteredState(1:indx,1)-obj.filteredState(1,1);
            y_all = obj.filteredState(1:indx,3);
            v_all = obj.filteredState(1:indx,6);
            a_all = obj.filteredState(1:indx,9);
            r = v_all./y_all;
            diffr = diff(r);
            
%             figure;
%             plot(y_all,r,'o');
%             text(y_all,r,string(1:length(r)));
%             title('y vs r')
            
            for ct=1:length(obj.rrefSegments)
                obj.rrefEntrySegments(ct) = data4rrefEntry();
                
                % Saving some parameters
                obj.rrefEntrySegments(ct).factor = obj.rrefSegments(ct).factor;
                
                % Find entry segment
                if isempty(obj.rrefSegments(ct).intervals_ti)
                    continue
                end
                
                
                rref_intervals = [1 1; sortrows(obj.rrefSegments(ct).intervals_ti,1,'ascend')];
                indices = nan(size(rref_intervals,1)-1,1);
                for ct1=2:size(rref_intervals,1)
                    if obj.rrefSegments(ct).rmean_ti(ct1-1) > -0.5
                        continue
                    end
                    
                    if diffr(rref_intervals(ct1,1)) > 0
                        indx1 = find(diffr(rref_intervals(ct1-1,2):rref_intervals(ct1,1)) < 0, 1, 'last');
                    elseif diffr(rref_intervals(ct1,1)) < 0
                        indx1 = find(diffr(rref_intervals(ct1-1,2):rref_intervals(ct1,1)) > 0, 1, 'last');
                    end
                    if isempty(indx1)
                        indx1 = rref_intervals(ct1-1,2);
                    else
                        indx1 = indx1 + rref_intervals(ct1-1,2);
                    end
                    
                    indx2 = find(r(rref_intervals(ct1-1,2):rref_intervals(ct1,1)) > -0.5, 1, 'last');
                    if isempty(indx2)
                        indx2 = rref_intervals(ct1-1,2);
                    else
                        indx2 = indx2 + rref_intervals(ct1-1,2);
                    end
                    
                    indx = max([indx1, indx2]);
                    
                    if rref_intervals(ct1,1)-indx>=15 && (obj.rrefSegments(ct).rmean_ti(ct1-1)*1.1 >= mean(r(indx:rref_intervals(ct1,1))) || ...
                                                          obj.rrefSegments(ct).rmean_ti(ct1-1)*0.9 <= mean(r(indx:rref_intervals(ct1,1))))
                        indices(ct1-1) = indx;
                    end
                end
                
                if size(rref_intervals,1)>1
                    [rref_intervals, sort_indices] = sortrows(obj.rrefSegments(ct).intervals_ti,1,'ascend');
                    nonnan_sort_indices = sort_indices(~isnan(indices));
                    obj.rrefEntrySegments(ct).intervals = [indices(~isnan(indices)) rref_intervals(~isnan(indices),:)]; % N by 2, where N is the number of non-overlapping time-window independent intervals
                    obj.rrefEntrySegments(ct).rref = obj.rrefSegments(ct).rref_ti(nonnan_sort_indices);
                    obj.rrefEntrySegments(ct).rmean = obj.rrefSegments(ct).rmean_ti(nonnan_sort_indices);
                end
            end
            
            
        end
        
        function find_rdot_estimate_in_rrefEntry(obj, start_end)
            % Goes through obj.rrefEntrySegments and finds rdot estimate
            % between [x1*rref x2*rref] entry segment where x1,
            % x2 vary from 0 to 1
            
            if isempty(obj.rrefEntrySegments)
                return
            end
            
            indx = size(obj.filteredState,1);
            t_all = obj.filteredState(1:indx,1)-obj.filteredState(end,1);
            y_all = obj.filteredState(1:indx,3);
            v_all = obj.filteredState(1:indx,6);
            a_all = obj.filteredState(1:indx,9);
            r = v_all./y_all;
%             diffr = diff(r);
            
%             figure;
%             plot(y_all,r,'o');
%             text(y_all,r,string(1:length(r)));
%             title('y vs r')
            
            for ct=1:length(obj.rrefEntrySegments)
                
                % Find entry segment
                rrefEntry_intervals = obj.rrefEntrySegments(ct).intervals;
                if isempty(rrefEntry_intervals)
                    continue
                end
                
                for ct1=1:size(rrefEntry_intervals,1)
                    rrefEntry_interval = rrefEntry_intervals(ct1,:);
                    rref = obj.rrefEntrySegments(ct).rref(ct1);
                    t_interval = t_all(rrefEntry_interval(1):rrefEntry_interval(2));
                    r_interval = r(rrefEntry_interval(1):rrefEntry_interval(2));   % it is negative (as y definition is negatie up until here)                 
                    y_interval = y_all(rrefEntry_interval(1):rrefEntry_interval(2));   % it is negative (as y definition is negatie up until here)                 
                    v_interval = v_all(rrefEntry_interval(1):rrefEntry_interval(2));   % it is negative (as y definition is negatie up until here)                 
                    a_interval = a_all(rrefEntry_interval(1):rrefEntry_interval(2));   % it is negative (as y definition is negatie up until here)                 
                    intercept_slope = [ones(length(t_interval),1) t_interval]\r_interval;
                    obj.rrefEntrySegments(ct).const_rvst(ct1,1) = intercept_slope(1);
                    obj.rrefEntrySegments(ct).slope_rvst(ct1,1) = intercept_slope(2);
                    
                    % Compute other parameters used for validation (R2,
                    % y_dist etc.)
                    rresid = r_interval - [ones(length(t_interval),1) t_interval]*intercept_slope;
                    SSresid = sum(rresid.^2); SStotal = (length(r_interval)-1)*var(r_interval);
                    obj.rrefEntrySegments(ct).Rsquare_rvst(ct1,1) = 1-SSresid/SStotal;
                    
                    obj.rrefEntrySegments(ct).ymean_for_rdot(ct1,1) = mean(y_all((rrefEntry_interval(1):rrefEntry_interval(2))));
                    obj.rrefEntrySegments(ct).delta_r(ct1,1) = abs(diff(r_interval([1 end])));
                    obj.rrefEntrySegments(ct).delta_y_actual(ct1,1) = abs(diff(y_interval([1 end])));
                    
                    tspan = [0 abs(diff(t_interval([1 end])))];
                    y0 = [-y_interval(1) v_interval(1)]; rdot = -intercept_slope(2);
                    [t,y] = ode45(@(t,y) filteredState_BlindLandingtrack.const_drdt_movement(t,y,rdot), tspan, y0);
                    obj.rrefEntrySegments(ct).delta_y_analytical(ct1,1) = abs(diff(y([1 end],1)));
                    asim = (rdot*y(:,1).^2-y(:,2).^2)./y(:,1);
                    
%                     y0 = [y_interval(1) v_interval(1)]; rdot = intercept_slope(2);
%                     [t,y] = ode45(@(t,y) filteredState_BlindLandingtrack.const_drdt_movement(t,y,rdot), tspan, y0);
%                     obj.rrefEntrySegments(ct).delta_y_analytical(ct1,1) = abs(diff(y([1 end],1)));
                    
%                     abs(diff(y([1 end],1))) - abs(diff(y_interval([1 end])))
                    
                    if -r_interval(1) < -r_interval(end)
                        obj.rrefEntrySegments(ct).isRise(ct1,1) = true;
                    elseif -r_interval(1) > -r_interval(end)
                        obj.rrefEntrySegments(ct).isRise(ct1,1) = false;
                    end
                    
                    obj.rrefEntrySegments(ct).yEntryStart(ct1,1) = y_interval(1);
                    obj.rrefEntrySegments(ct).vEntryStart(ct1,1) = v_interval(1);
                    obj.rrefEntrySegments(ct).rEntryStart(ct1,1) = r_interval(1);
                    obj.rrefEntrySegments(ct).delta_Ventry(ct1,1) = diff(v_interval([1 end]));
                    obj.rrefEntrySegments(ct).delta_tentry(ct1,1) = diff(t_interval([1 end]));
                    obj.rrefEntrySegments(ct).amean_entry(ct1,1) = mean(a_interval);
                end
            end
            
        end
        
        function [delta_t, dist_travelled] = compute_landing_performance(obj,y1,y2)
            % This function computes time taken by bbee to go from y1 to y2
            % It also computes trajectory path lengrh between y1 to y2
            
            assert(y2<y1);
            
            t = obj.filteredState(:,1)-obj.filteredState(1,1);
            y = -obj.filteredState(:,3);
            
            delta_t = nan;
            % Computing delta_t
            if y(1)>=y1
                indx_y1 = find((y-y1)<0, 1);
                indx_y2 = find((y-y2)<0, 1);
                
%                 if isempty(indx_y2)
%                     indx_y2 = find(abs(y-y2)<1e-3, 1);
%                 end                

                if isempty(indx_y1) || isempty(indx_y2)
                    keyboard;
                end
                assert(~isempty(indx_y1) & ~isempty(indx_y2));
                assert(indx_y2>indx_y1);
                delta_t = t(indx_y2)-t(indx_y1);                
            end
            
            % Distance travelled is not implemented yet
            dist_travelled = nan; 
        end
        
        function compute_instabilityFollows(obj, yrange)
            % For every rref segment, checks whether instability follows
            % within yrange after it ends. (Method 1)
            % Condition for instability (V<0.05 m/s)
            
            vthreshold = 0.05; % in m/s
            
            if isempty(obj.rrefSegments)
                return
            end
            
            indx = size(obj.filteredState,1);
            t_all = obj.filteredState(1:indx,1)-obj.filteredState(1,1);
            y_all = -obj.filteredState(1:indx,3);
            v_all = obj.filteredState(1:indx,6);
%             a_all = obj.filteredState(1:indx,9);
%             r = v_all./y_all;
            
            for ct=1:length(obj.rrefSegments)
                if isempty(obj.rrefSegments(ct).intervals_ti)
                    continue
                end
                rref_intervals = obj.rrefSegments(ct).intervals_ti;
                
%                 try
%                 dum = diff(rref_intervals);
%                 assert(all(dum(:)>0));
%                 catch
%                     keyboard
%                 end
%               Conclusion: rref segments are saved in increasing order of
%               indices i.e., time!
                
                instabilityFollows = false(size(rref_intervals,1),1);
                instability_indices = nan(size(rref_intervals,1),2);
                instability_deltat = nan(size(rref_intervals,1),1);
                instability_meanv = nan(size(rref_intervals,1),1);
                for ct1=1:size(rref_intervals,1)
                    
                    
                    y_part = y_all(rref_intervals(ct1,2):end)-y_all(rref_intervals(ct1,2));
                    indx = find(y_part<-yrange,1)+rref_intervals(ct1,2)-1; % where yrrefend+yrange occurs
                    
                    
                    if ct1<size(rref_intervals,1)
                        % search until the start of next rref segment or
                        % where yrrefend+yrange occurs (whichever occurs
                        % first)
                        indx = min([indx rref_intervals(ct1+1,1)]);  
                    elseif isempty(indx) && ct1==size(rref_intervals,1)
                        indx = min([indx length(y_all)]);
                    end
                                            
                    if any(v_all(rref_intervals(ct1,2):indx)<vthreshold)
                         
                        % After start of instability, find instant at which
                        % bbee reaches same y as at the start of
                        % instability with positive v. Then, look back and find instant at
                        % which velocity is smaller than vthreshold.
                        
                        % Where instability starts
                        indx1 = find(v_all(rref_intervals(ct1,2):indx)<vthreshold,1) ...
                                      +rref_intervals(ct1,2)-1; % 
                        
                        %
                        
                        indx2 = find(y_all(indx1:end)-y_all(indx1)<0 & ...
                            v_all(indx1:end)>vthreshold, 1)+indx1-1;
                        
                        indx2 = find(v_all(indx1:indx2)<vthreshold,1,'last')+indx1-1; % Find where instability ended
%                         [[1:length(y_all)]' y_all v_all]

%                         if y_all(indx2)>y_all(rref_intervals(ct1,1)) % instability ends further away from where rref starts
%                             continue; % bbee flew backwards
%                         try
                        if isempty(indx2) || y_all(indx2)>y_all(rref_intervals(ct1,1))
                            % First condition - bbee flew backwards,
                            % therefore instability end not stored
                            % Second condition - start found, not end 
                            instability_indices(ct1,1) = indx1;
                            instabilityFollows(ct1) = true;
                        else
%                             try
                            instability_indices(ct1,:) = [indx1 indx2];
                            instabilityFollows(ct1) = true;
                            instability_deltat(ct1) = diff(t_all([indx1 indx2]));
                            instability_meanv(ct1) = mean(v_all(indx1:indx2));
%                             catch
%                                 keyboard;
%                             end
                        end
%                         catch
%                             keyboard;
%                         end
%                         
                        
                    end
                    
                end
                obj.rrefSegments(ct).instabilityFollows = instabilityFollows;
                obj.rrefSegments(ct).instability_indices = instability_indices;
                obj.rrefSegments(ct).instability_deltat = instability_deltat;
                obj.rrefSegments(ct).instability_meanv = instability_meanv;
                obj.rrefSegments(ct).y_rrefEnd = arrayfun(@(x) y_all(x), rref_intervals(:,2));
            end
        end
        
        function compute_instabilityFollows2(obj)
            % For every rref segment, checks whether instability follows
            % i.e. if r decreases continuously after rref ends and V reaches < 0.05 m/s
            
            vthreshold = 0.05; % in m/s
            
            if isempty(obj.rrefSegments)
                return
            end
            
            indx = size(obj.filteredState,1);
            t_all = obj.filteredState(1:indx,1)-obj.filteredState(1,1);
            y_all = -obj.filteredState(1:indx,3);
            v_all = obj.filteredState(1:indx,6);
%             a_all = obj.filteredState(1:indx,9);
            r = v_all./y_all;
            diffr = diff(r);
            
            for ct=1:length(obj.rrefSegments)
                if isempty(obj.rrefSegments(ct).intervals_ti)
                    continue
                end
                rref_intervals = obj.rrefSegments(ct).intervals_ti;
                
%                 try
%                 dum = diff(rref_intervals);
%                 assert(all(dum(:)>0));
%                 catch
%                     keyboard
%                 end
%               Conclusion: rref segments are saved in increasing order of
%               indices i.e., time!
                
                instabilityFollows = false(size(rref_intervals,1),1);
                instability_indices = nan(size(rref_intervals,1),2);
                instability_deltat = nan(size(rref_intervals,1),1);
                instability_meanv = nan(size(rref_intervals,1),1);
                for ct1=1:size(rref_intervals,1)
                    
                    indx = find(diffr(rref_intervals(ct1,2):end) > 0, 1, 'first') + rref_intervals(ct1,2);
                    
%                     y_part = y_all(rref_intervals(ct1,2):end)-y_all(rref_intervals(ct1,2));
%                     indx = find(y_part<-yrange,1)+rref_intervals(ct1,2)-1; % where yrrefend+yrange occurs
%                     
%                     
%                     if ct1<size(rref_intervals,1)
%                         % search until the start of next rref segment or
%                         % where yrrefend+yrange occurs (whichever occurs
%                         % first)
%                         indx = min([indx rref_intervals(ct1+1,1)]);  
%                     elseif isempty(indx) && ct1==size(rref_intervals,1)
%                         indx = min([indx length(y_all)]);
%                     end
                                            
                    if any(v_all(rref_intervals(ct1,2):indx)<vthreshold)
                         
                        % After start of instability, find instant at which
                        % bbee reaches same y as at the start of
                        % instability with positive v. Then, look back and find instant at
                        % which velocity is smaller than vthreshold.
                        
                        % Where instability starts
                        indx1 = find(v_all(rref_intervals(ct1,2):indx)<vthreshold,1) ...
                                      +rref_intervals(ct1,2)-1; % 
                        
                        %
                        
                        indx2 = find(y_all(indx1:end)-y_all(indx1)<0 & ...
                            v_all(indx1:end)>vthreshold, 1)+indx1-1;
                        
                        indx2 = find(v_all(indx1:indx2)<vthreshold,1,'last')+indx1-1; % Find where instability ended
%                         [[1:length(y_all)]' y_all v_all]

%                         if y_all(indx2)>y_all(rref_intervals(ct1,1)) % instability ends further away from where rref starts
%                             continue; % bbee flew backwards
%                         try
                        if isempty(indx2) || y_all(indx2)>y_all(rref_intervals(ct1,1))
                            % First condition - bbee flew backwards,
                            % therefore instability end not stored
                            % Second condition - start found, not end 
                            instability_indices(ct1,1) = indx1;
                            instabilityFollows(ct1) = true;
                        else
%                             try
                            instability_indices(ct1,:) = [indx1 indx2];
                            instabilityFollows(ct1) = true;
                            instability_deltat(ct1) = diff(t_all([indx1 indx2]));
                            instability_meanv(ct1) = mean(v_all(indx1:indx2));
%                             catch
%                                 keyboard;
%                             end
                        end
%                         catch
%                             keyboard;
%                         end
%                         
                        
                    end
                    
                end
                obj.rrefSegments(ct).instabilityFollows = instabilityFollows;
                obj.rrefSegments(ct).instability_indices = instability_indices;
                obj.rrefSegments(ct).instability_deltat = instability_deltat;
                obj.rrefSegments(ct).instability_meanv = instability_meanv;
                obj.rrefSegments(ct).y_rrefEnd = arrayfun(@(x) y_all(x), rref_intervals(:,2));
            end
        end
        
        function output = hasLowV(obj, ygroups, vthreshold)
            % find if the current track has low V (< vthreshold) within
            % each ygroup
            % output - logical vector of size length(ygroup)-1 by 1
            assert(all(diff(ygroups)>0));
            
            y0 = 0.05; % neglecting the part when the bbee covers first 5 cm towards the platform
            y = -obj.filteredState(:,3);
            indx = find(y<y(1)-y0, 1, 'first');
            
            y = -obj.filteredState(indx:end,3);
            v = obj.filteredState(indx:end,6);
            
            output = arrayfun(@(x) any(v(y>=ygroups(x) & y<ygroups(x+1))<vthreshold), 1:length(ygroups)-1);
            
        end
        
        function plotHandles = plot_rrefsEntry(obj, factor)
            % Plots V vs y and V/y vs y highling change in r* within the same track
            % excerpt with positive V and positive y
            % This function plots time-window independent best intervals
            
            % factor - for which threshold factor plots are required to be
            % produced
            
            assert(~isempty(obj.rrefEntrySegments) && length(factor) == 1);
            ct_factor = find(abs([obj.rrefEntrySegments.factor] - factor) < 1e-6);
            if isempty(ct_factor)
                error('Can NOT find the r* entry intervals for the asked factor.');
            end
            intervals = obj.rrefEntrySegments(ct_factor).intervals;
            if isempty(intervals)
                plotHandles = [];
                return;
            end
            
            dataPerTrackExcerpt = struct.empty;
            for ct=1:size(intervals,1)
                dataPerTrackExcerpt(ct).state4rrefEstimate = obj.filteredState(intervals(ct,2):intervals(ct,3),:);
                dataPerTrackExcerpt(ct).state4rrefEntry = obj.filteredState(intervals(ct,1):intervals(ct,2),:);
                dataPerTrackExcerpt(ct).rref = obj.rrefEntrySegments(ct_factor).rref(ct);%mean(dataPerTrackExcerpt(ct).state4rrefEstimate(:,6)./dataPerTrackExcerpt(ct).state4rrefEstimate(:,3));
            end
            datacombined = [vertcat(dataPerTrackExcerpt.state4rrefEstimate); vertcat(dataPerTrackExcerpt.state4rrefEntry)];
            
            ymin = min(datacombined(:,3));
            [~, indx_ymin] = min(abs(obj.filteredState(:,3)-ymin));
            time_ymin = obj.filteredState(indx_ymin,1);
            % Start from 0.06s data earlier than dataTrackExcerpt.filteredState(indx_ymax,1)
            indx_start = find(obj.filteredState(:,1)-time_ymin<-0.06,1);
            if isempty(indx_start)
                indx_start = 1;
            end
            complete_state = obj.filteredState(indx_start:end,:);
            
            plotHandles = figure;
            p1 = subplot(2,1,1); hold on;
            plot(-complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('V (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3)], 'FontSize', 16);
            %     ylim([-0.2 p1.YLim(2)]);
            %     ylim([-0.2 0.6]);
            %     yticks([-0.2:0.2:0.6]);
            %     xticks([p1.XLim(1):0.1:0]);
            
            
            p2 = subplot(2,1,2); hold on;
            %     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('r (1/s)', 'FontSize', 16);
            %     ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
            xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            ylim([0 ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1))]);
            yticks(0:1:ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1)));
            %     xticks(p1.XLim(1):0.1:0);
            
            linkaxes([p1, p2],'x');
            if abs(floor((min(datacombined(:,3))/0.05))*0.05 - min(datacombined(:,3))) > 0.025
                xmin = floor((min(datacombined(:,3))/0.05))*0.05;
            else
                xmin = floor((min(datacombined(:,3))/0.05))*0.05-0.05;
            end
            xlim([0 -xmin]);
            %     xlim([0 Inf]);
            
            for ct=1:size(intervals,1)
                obj = dataPerTrackExcerpt(ct);
                state_subset = obj.state4rrefEstimate;
                state_4entry = obj.state4rrefEntry;
                figure(plotHandles(1))
                subplot(2,1,1);
                plot(-state_subset(:,3),state_subset(:,6),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(-state_4entry(:,3),state_4entry(:,6),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                plot(-state_subset(:,3),-obj.rref*-state_subset(:,3),'LineWidth',2,'Color',[69 117 180]./255);
                plot([0 -state_subset(end,3)],[0 -obj.rref*-state_subset(end,3)],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                subplot(2,1,2);
                plot(-state_subset(:,3),state_subset(:,6)./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(-state_4entry(:,3),state_4entry(:,6)./-state_4entry(:,3),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
                plot(-state_subset([1,end],3),[-obj.rref, -obj.rref],'LineWidth',2,'Color',[69 117 180]./255);
                plot([p1.XLim(1) -state_subset(1,3)],[-obj.rref -obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
            end
            
        end
        
        function plotHandles = plot_rrefsEntry_with_acc(obj, factor)
            % Plots V vs y and V/y vs y and a vs y highling change in r* within the same track
            % excerpt with positive V and positive y
            % This function plots time-window independent best intervals
            
            % factor - for which threshold factor plots are required to be
            % produced
            
            assert(~isempty(obj.rrefEntrySegments) && length(factor) == 1);
            ct_factor = find(abs([obj.rrefEntrySegments.factor] - factor) < 1e-6);
            if isempty(ct_factor)
                error('Can NOT find the r* entry intervals for the asked factor.');
            end
            intervals = obj.rrefEntrySegments(ct_factor).intervals;
            if isempty(intervals)
                plotHandles = [];
                return;
            end
            
            dataPerTrackExcerpt = struct.empty;
            for ct=1:size(intervals,1)
                dataPerTrackExcerpt(ct).state4rrefEstimate = obj.filteredState(intervals(ct,2):intervals(ct,3),:);
                dataPerTrackExcerpt(ct).state4rrefEntry = obj.filteredState(intervals(ct,1):intervals(ct,2),:);
                dataPerTrackExcerpt(ct).rref = obj.rrefEntrySegments(ct_factor).rref(ct);%mean(dataPerTrackExcerpt(ct).state4rrefEstimate(:,6)./dataPerTrackExcerpt(ct).state4rrefEstimate(:,3));
            end
            datacombined = [vertcat(dataPerTrackExcerpt.state4rrefEstimate); vertcat(dataPerTrackExcerpt.state4rrefEntry)];
            
            ymin = min(datacombined(:,3));
            [~, indx_ymin] = min(abs(obj.filteredState(:,3)-ymin));
            time_ymin = obj.filteredState(indx_ymin,1);
            % Start from 0.06s data earlier than dataTrackExcerpt.filteredState(indx_ymax,1)
            indx_start = find(obj.filteredState(:,1)-time_ymin<-0.06,1);
            if isempty(indx_start)
                indx_start = 1;
            end
            complete_state = obj.filteredState(indx_start:end,:);
            
            plotHandles = figure;
            p1 = subplot(3,1,1); hold on;
            plot(-complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('V (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3)], 'FontSize', 16);
            ylim([0 ceil(p1.YLim(2)*10)/10]);
            %     ylim([-0.2 0.6]);
            yticks([0:0.1:p1.YLim(2)]);
            %     xticks([p1.XLim(1):0.1:0]);
            
            
            p2 = subplot(3,1,2); hold on;
            %     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('r (1/s)', 'FontSize', 16);
            %     ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
%             xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            ylim([0 ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1))]);
            yticks(0:1:ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1)));
            %     xticks(p1.XLim(1):0.1:0);
            
            p3 = subplot(3,1,3); hold on;
            plot(-complete_state(:,3),complete_state(:,9),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('a (ms^{-2})', 'FontSize', 16);
            xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            yline(0,'--');
            ylim([-4 4]);
            yticks(-4:4:4);
            
            linkaxes([p1, p2, p3],'x');
            if abs(floor((min(datacombined(:,3))/0.05))*0.05 - min(datacombined(:,3))) > 0.025
                xmin = floor((min(datacombined(:,3))/0.05))*0.05;
            else
                xmin = floor((min(datacombined(:,3))/0.05))*0.05-0.05;
            end
            xlim([0 -xmin]);
            %     xlim([0 Inf]);
            
            for ct=1:size(intervals,1)
                obj = dataPerTrackExcerpt(ct);
                state_subset = obj.state4rrefEstimate;
                state_4entry = obj.state4rrefEntry;
                figure(plotHandles(1))
                subplot(3,1,1);
                plot(-state_subset(:,3),state_subset(:,6),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(-state_4entry(:,3),state_4entry(:,6),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                plot(-state_subset(:,3),-obj.rref*-state_subset(:,3),'LineWidth',2,'Color',[69 117 180]./255);
                plot([0 -state_subset(end,3)],[0 -obj.rref*-state_subset(end,3)],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                subplot(3,1,2);
                plot(-state_subset(:,3),state_subset(:,6)./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(-state_4entry(:,3),state_4entry(:,6)./-state_4entry(:,3),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
                plot(-state_subset([1,end],3),[-obj.rref, -obj.rref],'LineWidth',2,'Color',[69 117 180]./255);
                plot([p1.XLim(1) -state_subset(1,3)],[-obj.rref -obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                subplot(3,1,3);
                plot(-state_subset(:,3),state_subset(:,9),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(-state_4entry(:,3),state_4entry(:,9),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                
            end
            
        end
        
        function plotHandles = plot_rrefsEntry_with_rdotestimate(obj, factor)
            % Plots r vs t, r vs y and V vs y highling change in r* within the same track
            % excerpt with positive V and positive y
            % This function plots time-window independent best intervals
            
            % factor - for which threshold factor plots are required to be
            % produced
            
            assert(~isempty(obj.rrefEntrySegments) && length(factor) == 1);
            ct_factor = find(abs([obj.rrefEntrySegments.factor] - factor) < 1e-6);
            if isempty(ct_factor)
                error('Can NOT find the r* entry intervals for the asked factor.');
            end
            intervals = obj.rrefEntrySegments(ct_factor).intervals;
            if isempty(intervals)
                plotHandles = [];
                return;
            end
            
            dataPerTrackExcerpt = struct.empty;
            for ct=1:size(intervals,1)
                dataPerTrackExcerpt(ct).state4rrefEstimate = obj.filteredState(intervals(ct,2):intervals(ct,3),:);
                dataPerTrackExcerpt(ct).state4rrefEntry = obj.filteredState(intervals(ct,1):intervals(ct,2),:);
                dataPerTrackExcerpt(ct).rref = obj.rrefEntrySegments(ct_factor).rref(ct);%mean(dataPerTrackExcerpt(ct).state4rrefEstimate(:,6)./dataPerTrackExcerpt(ct).state4rrefEstimate(:,3));
                dataPerTrackExcerpt(ct).const_rvst = obj.rrefEntrySegments(ct_factor).const_rvst(ct);
                dataPerTrackExcerpt(ct).slope_rvst = obj.rrefEntrySegments(ct_factor).slope_rvst(ct);
            end
            datacombined = [vertcat(dataPerTrackExcerpt.state4rrefEstimate); vertcat(dataPerTrackExcerpt.state4rrefEntry)];
            
            ymin = min(datacombined(:,3));
            [~, indx_ymin] = min(abs(obj.filteredState(:,3)-ymin));
            time_ymin = obj.filteredState(indx_ymin,1);
            % Start from 0.06s data earlier than dataTrackExcerpt.filteredState(indx_ymax,1)
            indx_start = find(obj.filteredState(:,1)-time_ymin<-0.06,1);
            if isempty(indx_start)
                indx_start = 1;
            end
            complete_state = obj.filteredState(indx_start:end,:);
            
            plotHandles = figure;
            p0 = subplot(3,1,1); hold on;
            plot(complete_state(:,1)-complete_state(end,1),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('r (s-1)', 'FontSize', 16);
            xlabel('t (s)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on
            title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3)], 'FontSize', 16);
            ylim([0 ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1))]);
            yticks(0:1:ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1)));
            
            
            p1 = subplot(3,1,3); hold on;
            plot(-complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('V (m/s)', 'FontSize', 16);
            xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
%             title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3)], 'FontSize', 16);
            ylim([0 ceil(p1.YLim(2)*10)/10]);
            %     ylim([-0.2 0.6]);
            yticks([0:0.1:p1.YLim(2)]);
            %     xticks([p1.XLim(1):0.1:0]);
            
            
            p2 = subplot(3,1,2); hold on;
            %     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('r (s-1)', 'FontSize', 16);
            
            set(gca, 'FontSize', 16); %grid on;
            ylim([0 ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1))]);
            yticks(0:1:ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1)));
            %     xticks(p1.XLim(1):0.1:0);
            
            
            subplot(3,1,3);
            linkaxes([p1, p2],'x');
            if abs(floor((min(datacombined(:,3))/0.05))*0.05 - min(datacombined(:,3))) > 0.025
                xmin = floor((min(datacombined(:,3))/0.05))*0.05;
            else
                xmin = floor((min(datacombined(:,3))/0.05))*0.05-0.05;
            end
            xlim([0 -xmin]);
            %     xlim([0 Inf]);
            subplot(3,1,1);
%             indx = find(abs(complete_state(:,3)-xmin) < 5e-3, 1, 'last');
            ttt = datacombined(:,1)-complete_state(end,1);
%             xlim([floor(min(ttt)*10)/10 ceil(max(ttt)*10)/10]);
            xlim([floor(min(ttt)/0.5)*0.5 ceil(max(ttt)/0.5)*0.5]);
            xticks([floor(min(ttt)/0.5)*0.5:0.25:ceil(max(ttt)/0.5)*0.5]);
%             xlim([floor(min(ttt)*10)/10 0]);
%             xlim([floor(ttt(indx)*10)/10 0]);

            
            
            for ct=1:size(intervals,1)
                obj = dataPerTrackExcerpt(ct);
                state_subset = obj.state4rrefEstimate;
                state_4entry = obj.state4rrefEntry;
                figure(plotHandles(1))
                subplot(3,1,1);
                plot(state_subset(:,1)-complete_state(end,1),state_subset(:,6)./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(state_4entry(:,1)-complete_state(end,1),state_4entry(:,6)./-state_4entry(:,3),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                plot(state_subset([1,end],1)-complete_state(end,1),[-obj.rref, -obj.rref],'LineWidth',2,'Color',[69 117 180]./255);
%                 plot([p1.XLim(1) state_subset(1,1)-complete_state(end,1)],[-obj.rref -obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
                % Plot rdot estimate
                restimated = obj.const_rvst + obj.slope_rvst*(state_4entry(:,1)-complete_state(end,1));
                plot(state_4entry(:,1)-complete_state(end,1),-restimated,'LineWidth',2,'Color',[0 0 0]./255);

                
                
                subplot(3,1,3);
                plot(-state_subset(:,3),state_subset(:,6),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(-state_4entry(:,3),state_4entry(:,6),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                plot(-state_subset(:,3),-obj.rref*-state_subset(:,3),'LineWidth',2,'Color',[69 117 180]./255);
%                 plot([0 -state_subset(end,3)],[0 -obj.rref*-state_subset(end,3)],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                subplot(3,1,2);
                plot(-state_subset(:,3),state_subset(:,6)./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(-state_4entry(:,3),state_4entry(:,6)./-state_4entry(:,3),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
                plot(-state_subset([1,end],3),[-obj.rref, -obj.rref],'LineWidth',2,'Color',[69 117 180]./255);
%                 plot([p1.XLim(1) -state_subset(1,3)],[-obj.rref -obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                
            end
            
        end
        
        function plotHandles = plot_rdotSimulation_with_actualdata(obj, factor)
            % Plots y,v,a,r vs t for rdot simulation and actual data
            % factor - for which threshold factor plots are required to be
            % produced
            
            assert(~isempty(obj.rrefEntrySegments) && length(factor) == 1);
            ct_factor = find(abs([obj.rrefEntrySegments.factor] - factor) < 1e-6);
            if isempty(ct_factor)
                error('Can NOT find the r* entry intervals for the asked factor.');
            end
            rrefEntry_intervals = obj.rrefEntrySegments(ct_factor).intervals;
            if isempty(rrefEntry_intervals)
                plotHandles = [];
                return;
            end
            
            indx = size(obj.filteredState,1);
            t_all = obj.filteredState(1:indx,1)-obj.filteredState(end,1);
            y_all = obj.filteredState(1:indx,3);
            v_all = obj.filteredState(1:indx,6);
            a_all = obj.filteredState(1:indx,9);
            r = v_all./y_all;
            
            plotHandles = figure;
            p1 = subplot(4,1,1); hold on;
            p2 = subplot(4,1,2); hold on;
            p3 = subplot(4,1,3); hold on;
            p4 = subplot(4,1,4); hold on;
            
            for ct1=1:size(rrefEntry_intervals,1)
                rrefEntry_interval = rrefEntry_intervals(ct1,:);
                rref = obj.rrefEntrySegments(ct_factor).rref(ct1);
                t_interval = t_all(rrefEntry_interval(1):rrefEntry_interval(2));
                r_interval = r(rrefEntry_interval(1):rrefEntry_interval(2));   % it is negative (as y definition is negatie up until here)
                y_interval = y_all(rrefEntry_interval(1):rrefEntry_interval(2));   % it is negative (as y definition is negatie up until here)
                v_interval = v_all(rrefEntry_interval(1):rrefEntry_interval(2));   % it is negative (as y definition is negatie up until here)
                a_interval = a_all(rrefEntry_interval(1):rrefEntry_interval(2));   % it is negative (as y definition is negatie up until here)
                intercept_slope = [ones(length(t_interval),1) t_interval]\r_interval;
                
                
                tspan = [0 abs(diff(t_interval([1 end])))];
                y0 = [-y_interval(1) v_interval(1)]; rdot = -intercept_slope(2);
                [t,y] = ode45(@(t,y) filteredState_BlindLandingtrack.const_drdt_movement(t,y,rdot), tspan, y0);
                asim = (rdot*y(:,1).^2-y(:,2).^2)./y(:,1);
                figure(plotHandles)
                subplot(4,1,1);
                plot(t_interval,-y_interval,'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(t+t_interval(1),y(:,1),'LineWidth',2,'Color',[69 117 180]./255);
                subplot(4,1,2);
                plot(t_interval,v_interval,'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(t+t_interval(1),y(:,2),'LineWidth',2,'Color',[69 117 180]./255);
                subplot(4,1,3);
                plot(t_interval,-r_interval,'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(t+t_interval(1),y(:,2)./y(:,1),'LineWidth',2,'Color',[69 117 180]./255);
                subplot(4,1,4);
                plot(t_interval,a_interval,'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(t+t_interval(1),asim,'LineWidth',2,'Color',[69 117 180]./255);
            end
            
            figure(plotHandles)
            subplot(4,1,1); 
            ylabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16);
            subplot(4,1,2); 
            ylabel('V (ms-1)', 'FontSize', 16);
            set(gca, 'FontSize', 16);
            subplot(4,1,3);
            ylabel('r (s-1)', 'FontSize', 16);
            set(gca, 'FontSize', 16);
            subplot(4,1,4);
            ylabel('A (ms-2)', 'FontSize', 16);
            set(gca, 'FontSize', 16);
            
        end
        
        function plotHandles = plot_rrefsEntry_withSimulatedData(obj, factor, tf)
            % Plots V vs y and V/y vs y highling change in r* within the same track
            % excerpt with positive V and positive y
            % This function plots time-window independent best intervals
            
            % factor - chosen factor for which transfer function is
            % simulated
            
            % tf - model for which simulated data is to be plotted
            
            r = -obj.filteredState(:,6)./obj.filteredState(:,3);
            y = -obj.filteredState(:,3);
            
            assert(~isempty(obj.rrefEntrySegments) && length(factor) == 1);
            ct_factor = find(abs([obj.rrefEntrySegments.factor] - factor) < 1e-6);
            if isempty(ct_factor)
                error('Can NOT find the r* entry intervals for the asked factor.');
            end
            intervals = obj.rrefEntrySegments(ct_factor).intervals;
            if isempty(intervals)
                plotHandles = [];
                return;
            end
            
            dataPerTrackExcerpt = struct.empty;
            for ct=1:size(intervals,1)
                dataPerTrackExcerpt(ct).state4rrefEstimate = obj.filteredState(intervals(ct,2):intervals(ct,3),:);
                dataPerTrackExcerpt(ct).state4rrefEntry = obj.filteredState(intervals(ct,1):intervals(ct,2),:);
                dataPerTrackExcerpt(ct).rref = obj.rrefEntrySegments(ct_factor).rref(ct);%mean(dataPerTrackExcerpt(ct).state4rrefEstimate(:,6)./dataPerTrackExcerpt(ct).state4rrefEstimate(:,3));
            end
            datacombined = [vertcat(dataPerTrackExcerpt.state4rrefEstimate); vertcat(dataPerTrackExcerpt.state4rrefEntry)];
            
            ymin = min(datacombined(:,3));
            [~, indx_ymin] = min(abs(obj.filteredState(:,3)-ymin));
            time_ymin = obj.filteredState(indx_ymin,1);
            % Start from 0.06s data earlier than dataTrackExcerpt.filteredState(indx_ymax,1)
            indx_start = find(obj.filteredState(:,1)-time_ymin<-0.06,1);
            if isempty(indx_start)
                indx_start = 1;
            end
            complete_state = obj.filteredState(indx_start:end,:);
            
            plotHandles = figure;
            p1 = subplot(2,1,1); hold on;
            plot(-complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('V (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            %     ylim([-0.2 p1.YLim(2)]);
            %     ylim([-0.2 0.6]);
            %     yticks([-0.2:0.2:0.6]);
            %     xticks([p1.XLim(1):0.1:0]);
            
            
            p2 = subplot(2,1,2); hold on;
            %     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('r (1/s)', 'FontSize', 16);
            %     ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
            xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            ylim([0 ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1))]);
            yticks(0:1:ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1)));
            %     xticks(p1.XLim(1):0.1:0);
            
            linkaxes([p1, p2],'x');
            if abs(floor((min(datacombined(:,3))/0.05))*0.05 - min(datacombined(:,3))) > 0.025
                xmin = floor((min(datacombined(:,3))/0.05))*0.05;
            else
                xmin = floor((min(datacombined(:,3))/0.05))*0.05-0.05;
            end
            xlim([0 -xmin]);
            %     xlim([0 Inf]);
            
            fitPercents = [];
            for ct=1:size(intervals,1)
                obj = dataPerTrackExcerpt(ct);
                state_subset = obj.state4rrefEstimate;
                state_4entry = obj.state4rrefEntry;
                
                output = r(intervals(ct,1):intervals(ct,3));
                input = -obj.rref*ones(length(output),1);
                dt = mean(diff(state_subset(:,1)));
                data1 = iddata(output, input, dt);

                if length(tf)==1
                    [rsim,fit,rsim0] = compare(data1, tf);
                else
                    assert(length(tf) == size(intervals,1));
                    [rsim,fit,rsim0] = compare(data1, tf(:,:,ct));
                end
                
                figure(plotHandles(1))
                subplot(2,1,1);
                plot(-state_subset(:,3),state_subset(:,6),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(-state_4entry(:,3),state_4entry(:,6),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                plot(y(intervals(ct,1):intervals(ct,3)),rsim.y.*y(intervals(ct,1):intervals(ct,3)),'LineWidth',2,'Color',[69 117 180]./255);
                plot([0 -state_subset(end,3)],[0 -obj.rref*-state_subset(end,3)],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                subplot(2,1,2);
                plot(-state_subset(:,3),state_subset(:,6)./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(-state_4entry(:,3),state_4entry(:,6)./-state_4entry(:,3),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
                plot(y(intervals(ct,1):intervals(ct,3)),rsim.y,'LineWidth',2,'Color',[69 117 180]./255);
                plot([p1.XLim(1) -state_subset(1,3)],[-obj.rref -obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                fitPercents = [fitPercents fit];
            end
            figure(plotHandles(1))
            subplot(2,1,1);
            title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3) ', % fit: ' num2str(fitPercents,4)], 'FontSize', 16);

            
        end
        
        function plotHandles = plot_rrefsEntry_actual_vs_filtered(obj, factor, fc)
            % Plots V vs y and V/y vs y highling entry in r* and r* within the same track
            % excerpt with positive V and positive y
            % This function plots time-window independent best intervals
            
            % factor - chosen factor for which plot is desired
            
            % fc - cut-off frequency for prefiltering using butterworth filter
            
            r = -obj.filteredState(:,6)./obj.filteredState(:,3);
            y = -obj.filteredState(:,3);
            
            assert(~isempty(obj.rrefEntrySegments) && length(factor) == 1);
            ct_factor = find(abs([obj.rrefEntrySegments.factor] - factor) < 1e-6);
            if isempty(ct_factor)
                error('Can NOT find the r* entry intervals for the asked factor.');
            end
            intervals = obj.rrefEntrySegments(ct_factor).intervals;
            if isempty(intervals)
                plotHandles = [];
                return;
            end
            
            dataPerTrackExcerpt = struct.empty;
            for ct=1:size(intervals,1)
                dataPerTrackExcerpt(ct).state4rrefEstimate = obj.filteredState(intervals(ct,2):intervals(ct,3),:);
                dataPerTrackExcerpt(ct).state4rrefEntry = obj.filteredState(intervals(ct,1):intervals(ct,2),:);
                dataPerTrackExcerpt(ct).stateEntryRrefCombined = obj.filteredState(intervals(ct,1):intervals(ct,3),:);
                dataPerTrackExcerpt(ct).rref = obj.rrefEntrySegments(ct_factor).rref(ct);%mean(dataPerTrackExcerpt(ct).state4rrefEstimate(:,6)./dataPerTrackExcerpt(ct).state4rrefEstimate(:,3));
            end
            datacombined = [vertcat(dataPerTrackExcerpt.state4rrefEstimate); vertcat(dataPerTrackExcerpt.state4rrefEntry)];
            
            ymin = min(datacombined(:,3));
            [~, indx_ymin] = min(abs(obj.filteredState(:,3)-ymin));
            time_ymin = obj.filteredState(indx_ymin,1);
            % Start from 0.06s data earlier than dataTrackExcerpt.filteredState(indx_ymax,1)
            indx_start = find(obj.filteredState(:,1)-time_ymin<-0.06,1);
            if isempty(indx_start)
                indx_start = 1;
            end
            complete_state = obj.filteredState(indx_start:end,:);
            
            plotHandles = figure;
            p1 = subplot(2,1,1); hold on;
            plot(-complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('V (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            %     ylim([-0.2 p1.YLim(2)]);
            %     ylim([-0.2 0.6]);
            %     yticks([-0.2:0.2:0.6]);
            %     xticks([p1.XLim(1):0.1:0]);
            
            
            p2 = subplot(2,1,2); hold on;
            %     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('r (1/s)', 'FontSize', 16);
            %     ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
            xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            ylim([0 ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1))]);
            yticks(0:1:ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1)));
            %     xticks(p1.XLim(1):0.1:0);
            
            linkaxes([p1, p2],'x');
            if abs(floor((min(datacombined(:,3))/0.05))*0.05 - min(datacombined(:,3))) > 0.025
                xmin = floor((min(datacombined(:,3))/0.05))*0.05;
            else
                xmin = floor((min(datacombined(:,3))/0.05))*0.05-0.05;
            end
            xlim([0 -xmin]);
            %     xlim([0 Inf]);
            
%             output = -state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),6)./...
%                 state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),3);
%             input = -rrefEntrySegments_fac.rmean(ct2)*ones(length(output),1);
%             dt = mean(diff(state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),1)));
                        
            fitPercents = [];
            for ct=1:size(intervals,1)
                obj = dataPerTrackExcerpt(ct);
                state_subset = obj.state4rrefEstimate;
                state_4entry = obj.state4rrefEntry;
                state_entryRrefCombined = obj.stateEntryRrefCombined;
                
                output = r(intervals(ct,1):intervals(ct,3));
                input = -obj.rref*ones(length(output),1);
                dt = mean(diff(state_subset(:,1)));
%                 data1 = iddata(output, input, dt);
                
                [num, den] = butter(2, fc/(1/dt/2),'low');
                rFiltered = filtfilt(num, den, output);
                
%                 if length(tf)==1
%                     [rsim,fit,rsim0] = compare(data1, tf);
%                 else
%                     assert(length(tf) == size(intervals,1));
%                     [rsim,fit,rsim0] = compare(data1, tf(:,:,ct));
%                 end
                
                figure(plotHandles(1))
                subplot(2,1,1);
                plot(-state_subset(:,3),state_subset(:,6),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(-state_4entry(:,3),state_4entry(:,6),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
%                 plot(y(intervals(ct,1):intervals(ct,3)),rsim.y.*y(intervals(ct,1):intervals(ct,3)),'LineWidth',2,'Color',[69 117 180]./255);
                plot(y(intervals(ct,1):intervals(ct,3)),rFiltered.*y(intervals(ct,1):intervals(ct,3)),'LineWidth',2,'Color',[169 169 169]./255);
                plot([0 -state_subset(end,3)],[0 -obj.rref*-state_subset(end,3)],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                subplot(2,1,2);
                plot(-state_subset(:,3),state_subset(:,6)./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(-state_4entry(:,3),state_4entry(:,6)./-state_4entry(:,3),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
%                 plot(y(intervals(ct,1):intervals(ct,3)),rsim.y,'LineWidth',2,'Color',[69 117 180]./255);
                plot(y(intervals(ct,1):intervals(ct,3)),rFiltered,'LineWidth',2,'Color',[169 169 169]./255);
                plot([p1.XLim(1) -state_subset(1,3)],[-obj.rref -obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
                
%                 fitPercents = [fitPercents fit];
            end
            figure(plotHandles(1))
            subplot(2,1,1);
%             title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3) ', % fit: ' num2str(fitPercents,4)], 'FontSize', 16);

            
        end
        
        function plotHandles = plot_rrefsEntry_withActualFilteredEstimatedData(obj, factor, sysiddata, tf)
            % Plots V vs y and V/y vs y highling entry in r* and r* within the same track
            % excerpt with positive V and positive y
            % This function plots time-window independent best intervals
            
            % factor - chosen factor for which plot is desired
            
            % sysiddata - data used for sys id (in cell array)
            % tf - transfer function array (3D array)
            
            r = -obj.filteredState(:,6)./obj.filteredState(:,3);
            y = -obj.filteredState(:,3);
            
            assert(~isempty(obj.rrefEntrySegments) && length(factor) == 1);
            
            ct_factor = find(abs([obj.rrefEntrySegments.factor] - factor) < 1e-6);
            if isempty(ct_factor)
                error('Can NOT find the r* entry intervals for the asked factor.');
            end
            intervals = obj.rrefEntrySegments(ct_factor).intervals;
            if isempty(intervals)
                plotHandles = [];
                return;
            end
            assert(size(intervals,1)==length(sysiddata) )%&& size(intervals,1)==size(tfs,3));
            
            dataPerTrackExcerpt = struct.empty;
            for ct=1:size(intervals,1)
                dataPerTrackExcerpt(ct).state4rrefEstimate = obj.filteredState(intervals(ct,2):intervals(ct,3),:);
                dataPerTrackExcerpt(ct).state4rrefEntry = obj.filteredState(intervals(ct,1):intervals(ct,2),:);
                dataPerTrackExcerpt(ct).stateEntryRrefCombined = obj.filteredState(intervals(ct,1):intervals(ct,3),:);
                dataPerTrackExcerpt(ct).rref = obj.rrefEntrySegments(ct_factor).rref(ct);%mean(dataPerTrackExcerpt(ct).state4rrefEstimate(:,6)./dataPerTrackExcerpt(ct).state4rrefEstimate(:,3));
            end
            datacombined = [vertcat(dataPerTrackExcerpt.state4rrefEstimate); vertcat(dataPerTrackExcerpt.state4rrefEntry)];
            
            ymin = min(datacombined(:,3));
            [~, indx_ymin] = min(abs(obj.filteredState(:,3)-ymin));
            time_ymin = obj.filteredState(indx_ymin,1);
            % Start from 0.06s data earlier than dataTrackExcerpt.filteredState(indx_ymax,1)
            indx_start = find(obj.filteredState(:,1)-time_ymin<-0.06,1);
            if isempty(indx_start)
                indx_start = 1;
            end
            complete_state = obj.filteredState(indx_start:end,:);
            
            plotHandles = figure;
            p1 = subplot(2,1,1); hold on;
            plot(-complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('V (ms-1)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            %     ylim([-0.2 p1.YLim(2)]);
            %     ylim([-0.2 0.6]);
            %     yticks([-0.2:0.2:0.6]);
            %     xticks([p1.XLim(1):0.1:0]);
            
            
            p2 = subplot(2,1,2); hold on;
            %     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('r (s-1)', 'FontSize', 16);
            %     ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
            xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            ylim([0 ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1))]);
            yticks(0:1:ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1)));
            %     xticks(p1.XLim(1):0.1:0);
            
            linkaxes([p1, p2],'x');
            if abs(floor((min(datacombined(:,3))/0.05))*0.05 - min(datacombined(:,3))) > 0.025
                xmin = floor((min(datacombined(:,3))/0.05))*0.05;
            else
                xmin = floor((min(datacombined(:,3))/0.05))*0.05-0.05;
            end
            xlim([0 -xmin]);
            %     xlim([0 Inf]);
            
%             output = -state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),6)./...
%                 state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),3);
%             input = -rrefEntrySegments_fac.rmean(ct2)*ones(length(output),1);
%             dt = mean(diff(state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),1)));
                        
            fitPercents = [];
            for ct=1:size(intervals,1)
                obj = dataPerTrackExcerpt(ct);
                state_subset = obj.state4rrefEstimate;
                state_4entry = obj.state4rrefEntry;
                state_entryRrefCombined = obj.stateEntryRrefCombined;
                
%                 output = r(intervals(ct,1):intervals(ct,3));
%                 input = -obj.rref*ones(length(output),1);
%                 dt = mean(diff(state_subset(:,1)));
%                 data1 = iddata(output, input, dt);
                
%                 [num, den] = butter(2, fc/(1/dt/2),'low');
                rFiltered = sysiddata{ct}.y;
                
                if length(tf)==1
                    [rsim,fit,rsim0] = compare(sysiddata{ct}, tf);
                else
                    assert(length(tf) == size(intervals,1));
                    [rsim,fit,rsim0] = compare(sysiddata{ct}, tf(:,:,ct));
                end
                
                figure(plotHandles(1))
                subplot(2,1,1);
                plot(-state_4entry(:,3),state_4entry(:,6),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                plot(-state_subset(:,3),state_subset(:,6),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                
                plot(y(intervals(ct,1):intervals(ct,3)),rFiltered.*y(intervals(ct,1):intervals(ct,3)),'LineWidth',2.5,'Color',[255 0 255]./255);
                plot(y(intervals(ct,1):intervals(ct,3)),rsim.y.*y(intervals(ct,1):intervals(ct,3)),'LineWidth',1.5,'Color',[69 117 180]./255);
                
                plot([0 -state_subset(1,3)],[0 -obj.rref*-state_subset(1,3)],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                subplot(2,1,2);
                plot(-state_4entry(:,3),state_4entry(:,6)./-state_4entry(:,3),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                plot(-state_subset(:,3),state_subset(:,6)./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                
                %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
                
                plot(y(intervals(ct,1):intervals(ct,3)),rFiltered,'LineWidth',2.5,'Color',[255 0 255]./255);
                plot(y(intervals(ct,1):intervals(ct,3)),rsim.y,'LineWidth',1.5,'Color',[69 117 180]./255);
                
                plot([p1.XLim(1) -state_subset(1,3)],[-obj.rref -obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                fitPercents = [fitPercents fit];
            end
            figure(plotHandles(1))
            subplot(2,1,1);
            title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3) ', % fit: ' num2str(fitPercents,4)], 'FontSize', 16);

            
        end
        
        function plotHandles = plot_rrefsEntry_withActualFilteredEstimatedData2(obj, factor, sysiddata, tf)
            % Plots r vs t, r vs y and V vs y 
            % Also highlights entry in r* and r* within the same track
            % excerpt with positive V and positive y
            % Overlays low-pass filtered data and sys-id data as well
            % This function plots time-window independent best intervals
            
            % factor - chosen factor for which plot is desired
            
            % sysiddata - data used for sys id (in cell array)
            % tf - transfer function array (3D array)
            
            r = -obj.filteredState(:,6)./obj.filteredState(:,3);
            y = -obj.filteredState(:,3);
            t = obj.filteredState(:,1)-obj.filteredState(end,1);
            
            assert(~isempty(obj.rrefEntrySegments) && length(factor) == 1);
            
            ct_factor = find(abs([obj.rrefEntrySegments.factor] - factor) < 1e-6);
            if isempty(ct_factor)
                error('Can NOT find the r* entry intervals for the asked factor.');
            end
            intervals = obj.rrefEntrySegments(ct_factor).intervals;
            if isempty(intervals)
                plotHandles = [];
                return;
            end
            assert(size(intervals,1)==length(sysiddata) )%&& size(intervals,1)==size(tfs,3));
            
            dataPerTrackExcerpt = struct.empty;
            for ct=1:size(intervals,1)
                dataPerTrackExcerpt(ct).state4rrefEstimate = obj.filteredState(intervals(ct,2):intervals(ct,3),:);
                dataPerTrackExcerpt(ct).state4rrefEntry = obj.filteredState(intervals(ct,1):intervals(ct,2),:);
                dataPerTrackExcerpt(ct).stateEntryRrefCombined = obj.filteredState(intervals(ct,1):intervals(ct,3),:);
                dataPerTrackExcerpt(ct).rref = obj.rrefEntrySegments(ct_factor).rref(ct);%mean(dataPerTrackExcerpt(ct).state4rrefEstimate(:,6)./dataPerTrackExcerpt(ct).state4rrefEstimate(:,3));
            end
            datacombined = [vertcat(dataPerTrackExcerpt.state4rrefEstimate); vertcat(dataPerTrackExcerpt.state4rrefEntry)];
            
            ymin = min(datacombined(:,3));
            [~, indx_ymin] = min(abs(obj.filteredState(:,3)-ymin));
            time_ymin = obj.filteredState(indx_ymin,1);
            % Start from 0.06s data earlier than dataTrackExcerpt.filteredState(indx_ymax,1)
            indx_start = find(obj.filteredState(:,1)-time_ymin<-0.06,1);
            if isempty(indx_start)
                indx_start = 1;
            end
            complete_state = obj.filteredState(indx_start:end,:);
            
            plotHandles = figure;
            p0 = subplot(3,1,1); hold on;
            plot(complete_state(:,1)-complete_state(end,1),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('r (s-1)', 'FontSize', 16);
            xlabel('t (s)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on
            title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3)], 'FontSize', 16);
            ylim([0 ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1))]);
            yticks(0:1:ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1)));
            
            p1 = subplot(3,1,3); hold on;
            plot(-complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('V (ms-1)', 'FontSize', 16);
            xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            %     ylim([-0.2 p1.YLim(2)]);
            %     ylim([-0.2 0.6]);
            yticks([floor(p1.YLim(1)):0.1:ceil(p1.YLim(2))]);
            %     xticks([p1.XLim(1):0.1:0]);
            
            
            p2 = subplot(3,1,2); hold on;
            %     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('r (s-1)', 'FontSize', 16);
            %     ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            ylim([0 ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1))]);
            yticks(0:1:ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1)));
            %     xticks(p1.XLim(1):0.1:0);
            
            subplot(3,1,3)
            linkaxes([p1, p2],'x');
            if abs(floor((min(datacombined(:,3))/0.05))*0.05 - min(datacombined(:,3))) > 0.025
                xmin = floor((min(datacombined(:,3))/0.05))*0.05;
            else
                xmin = floor((min(datacombined(:,3))/0.05))*0.05-0.05;
            end
            xlim([0 -xmin]);
            %     xlim([0 Inf]);
            subplot(3,1,1);
%             indx = find(abs(complete_state(:,3)-xmin) < 5e-3, 1, 'last');
            ttt = datacombined(:,1)-complete_state(end,1);
%             xlim([floor(min(ttt)*10)/10 ceil(max(ttt)*10)/10]);
            xlim([floor(min(ttt)/0.5)*0.5 ceil(max(ttt)/0.5)*0.5]);
            xticks([floor(min(ttt)/0.5)*0.5:0.25:ceil(max(ttt)/0.5)*0.5]);
            
%             output = -state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),6)./...
%                 state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),3);
%             input = -rrefEntrySegments_fac.rmean(ct2)*ones(length(output),1);
%             dt = mean(diff(state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),1)));
                        
            fitPercents = [];
            for ct=1:size(intervals,1)
                obj = dataPerTrackExcerpt(ct);
                state_subset = obj.state4rrefEstimate;
                state_4entry = obj.state4rrefEntry;
                state_entryRrefCombined = obj.stateEntryRrefCombined;
                
%                 output = r(intervals(ct,1):intervals(ct,3));
%                 input = -obj.rref*ones(length(output),1);
%                 dt = mean(diff(state_subset(:,1)));
%                 data1 = iddata(output, input, dt);
                
%                 [num, den] = butter(2, fc/(1/dt/2),'low');
                rFiltered = sysiddata{ct}.y;
                
                if length(tf)==1
                    [rsim,fit,rsim0] = compare(sysiddata{ct}, tf);
                else
                    assert(length(tf) == size(intervals,1));
                    [rsim,fit,rsim0] = compare(sysiddata{ct}, tf(:,:,ct));
                end
                
                figure(plotHandles(1))
                subplot(3,1,1);
                plot(state_subset(:,1)-complete_state(end,1),state_subset(:,6)./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(state_4entry(:,1)-complete_state(end,1),state_4entry(:,6)./-state_4entry(:,3),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                plot(state_subset([1,end],1)-complete_state(end,1),[-obj.rref, -obj.rref],'LineWidth',2,'Color',[69 117 180]./255);
                plot(t(intervals(ct,1):intervals(ct,3)),rFiltered,'LineWidth',2.5,'Color',[255 0 255]./255);
                plot(t(intervals(ct,1):intervals(ct,3)),rsim.y,'LineWidth',1.5,'Color',[69 117 180]./255);
                
%                 plot([p1.XLim(1) -state_subset(1,3)],[-obj.rref -obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                subplot(3,1,3);
                plot(-state_4entry(:,3),state_4entry(:,6),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                plot(-state_subset(:,3),state_subset(:,6),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                
                plot(y(intervals(ct,1):intervals(ct,3)),rFiltered.*y(intervals(ct,1):intervals(ct,3)),'LineWidth',2.5,'Color',[255 0 255]./255);
                plot(y(intervals(ct,1):intervals(ct,3)),rsim.y.*y(intervals(ct,1):intervals(ct,3)),'LineWidth',1.5,'Color',[69 117 180]./255);
                
                plot([0 -state_subset(1,3)],[0 -obj.rref*-state_subset(1,3)],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                subplot(3,1,2);
                plot(-state_4entry(:,3),state_4entry(:,6)./-state_4entry(:,3),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
                plot(-state_subset(:,3),state_subset(:,6)./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                
                %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
                
                plot(y(intervals(ct,1):intervals(ct,3)),rFiltered,'LineWidth',2.5,'Color',[255 0 255]./255);
                plot(y(intervals(ct,1):intervals(ct,3)),rsim.y,'LineWidth',1.5,'Color',[69 117 180]./255);
                
                plot([p1.XLim(1) -state_subset(1,3)],[-obj.rref -obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                fitPercents = [fitPercents fit];
            end
            figure(plotHandles(1))
            subplot(3,1,1);
            title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3) ', % fit: ' num2str(fitPercents,4)], 'FontSize', 16);

            
        end
    
        
        function plotHandles = plot_rrefsEntry_with_time(obj, factor)
            % Plots V vs t and V/y vs t highling change in r* within the same track
            % excerpt with positive V
            % This function plots time-window independent best intervals
            
            % factor - for which threshold factor plots are required to be
            % produced
            
            assert(~isempty(obj.rrefEntrySegments) && length(factor) == 1);
            ct_factor = find(abs([obj.rrefEntrySegments.factor] - factor) < 1e-6);
            if isempty(ct_factor)
                error('Can NOT find the r* entry intervals for the asked factor.');
            end
            intervals = obj.rrefEntrySegments(ct_factor).intervals;
            if isempty(intervals)
                plotHandles = [];
                return;
            end
            
            dataPerTrackExcerpt = struct.empty;
            for ct=1:size(intervals,1)
                dataPerTrackExcerpt(ct).state4rrefEstimate = obj.filteredState(intervals(ct,2):intervals(ct,3),:);
                dataPerTrackExcerpt(ct).state4rrefEntry = obj.filteredState(intervals(ct,1):intervals(ct,2),:);
                dataPerTrackExcerpt(ct).rref = obj.rrefEntrySegments(ct_factor).rref(ct);%mean(dataPerTrackExcerpt(ct).state4rrefEstimate(:,6)./dataPerTrackExcerpt(ct).state4rrefEstimate(:,3));
            end
            datacombined = [vertcat(dataPerTrackExcerpt.state4rrefEstimate); vertcat(dataPerTrackExcerpt.state4rrefEntry)];
            
            ymin = min(datacombined(:,3));
            [~, indx_ymin] = min(abs(obj.filteredState(:,3)-ymin));
            time_ymin = obj.filteredState(indx_ymin,1);
            % Start from 0.06s data earlier than dataTrackExcerpt.filteredState(indx_ymax,1)
            indx_start = find(obj.filteredState(:,1)-time_ymin<-0.06,1);
            if isempty(indx_start)
                indx_start = 1;
            end
            complete_state = obj.filteredState(indx_start:end,:);
            tend = obj.filteredState(end,1);
            
            plotHandles = figure; hold on;
%             plot(-complete_state(:,1)-tend,complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('r (1/s)', 'FontSize', 16);
            %     ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
            xlabel('t (s)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            ylim([0 ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1))]);
            yticks(0:1:ceil(-1*ceil(min(datacombined(:,6)./datacombined(:,3))-1)));
            
%             if abs(floor((min(datacombined(:,3))/0.05))*0.05 - min(datacombined(:,3))) > 0.025
%                 xmin = floor((min(datacombined(:,3))/0.05))*0.05;
%             else
%                 xmin = floor((min(datacombined(:,3))/0.05))*0.05-0.05;
%             end
%             xlim([0 -xmin]);
            %     xlim([0 Inf]);
            
            for ct=1:size(intervals,1)
                obj = dataPerTrackExcerpt(ct);
                state_subset = obj.state4rrefEstimate;
                state_4entry = obj.state4rrefEntry;
                figure(plotHandles(1))
                plot(state_subset(:,1)-tend,state_subset(:,6)./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(state_4entry(:,1)-tend,state_4entry(:,6)./-state_4entry(:,3),'.','MarkerSize',12,'MarkerFaceColor',[161 217 155]./255, 'MarkerEdgeColor',[161 217 155]./255);
            end
            
        end

        
%         function compute_rref1(obj, params, factors, time_window)
%             % params - [sigma_{rmean-c_{r vs y, 0-1}}, sigma_{m_{r vs y, 0-1}}] = [sigma1 sigma2]
%             % factors - vector of some size
%             % time_window - 2 by 1 vector = [min_gap max_gap]
%             
%             % For each point, intervals are found by looking ahead
%             % min_gap:1:max_gap points
%             
%             % Within each such interval 3 lines are fit
%             % First line between whole interval: r = m_{0-1} y + c_{0-1}
%             % Second line between first half of the interval: r = m_{0-0.5} y + c_{0-0.5}
%             % Third line between second half of the interval: r = m_{0.5-1} y + c_{0.5-1}
%             % Time windows/intervals which satisfy the following constraint are
%             % classsified as intervals with constant r (or simply r*).
%             % Constraint: 
%             
%             assert(numel(params)==2 && numel(time_window)==2);
%             rmse_func = @(idx, indices, data) sqrt(sum((data(indices(idx,1):indices(idx,2)) - mean(data(indices(idx,1):indices(idx,2)))).^2) ...
%                                                    /(indices(idx,2)-indices(idx,1)+1));
%             
%             min_gap = time_window(1); % the first straight line is between current point and (current point + min_gap)th point
%             max_gap = time_window(2);
%             
%             N = size(obj.filteredState,1);
%             
%             t_all = obj.filteredState(:,1)-obj.filteredState(1,1);
%             y_all = obj.filteredState(:,3);
%             v_all = obj.filteredState(:,6);
%             a_all = obj.filteredState(:,9);
%             r_all = v_all./y_all;
% 
%             x_all = obj.filteredState(:,2);
%             z_all = obj.filteredState(:,4);
%             
% %             figure;
% %             plot(t_all,r_all,'o');
% %             text(t_all,r_all,string(1:length(r_all)));
% %             title('t vs r')
%             
% %             figure;
% %             plot(y_all,r_all,'o');
% %             text(y_all,r_all,string(1:length(r_all)));
% %             title('y vs r')
%             
% %             figure;
% %             plot(t_all,rdot_all,'ro');
% %             text(t_all,rdot_all,string(1:length(rdot_all)));
% %             title('t vs rdot');
%             
%             for ct_factor=1:length(factors)
%                 obj.rrefSegments(ct_factor) = data4rrefEstimate();
%                 
%                 % Saving some parameters
%                 obj.rrefSegments(ct_factor).factor = factors(ct_factor);
%                 obj.rrefSegments(ct_factor).params = params;
%             end
%             
%             possible_indices = cell(length(factors),1); % possible intervals for each factor
%             
%             for ct=1:N-min_gap % for each point as a starting point for different lines through origin
%                 if v_all(ct)/y_all(ct) >= 0
%                     continue;
%                 end
%                 
%                 end_point = min([ct+max_gap, N]);
%                 
%                 % Useful column vectors
% %                 t = t_all(ct:end_point)-t_all(1);
%                 y = y_all(ct:end_point);
%                 v = v_all(ct:end_point);
%                 a = a_all(ct:end_point);
%                 r = v./y;
% %                 rdot = abs(a./y-(v./y).^2);
%                 
%                 last_index = find(r > 0, 1) - 1; % last but one point where r becomes positive
%                 if isempty(last_index)
%                     last_index = length(r);
%                 end
%                 
%                 coeff = zeros(last_index-min_gap,2); % [intercept slope]
%                 coeff_firsthalf = zeros(last_index-min_gap,2); 
%                 coeff_secondhalf = zeros(last_index-min_gap,2); 
%                 r_mean = zeros(last_index-min_gap,1); % Mean of r within each interval
%                 a_mean = zeros(last_index-min_gap,1); % Mean of a within each interval
%                 a_mean_firsthalf = zeros(last_index-min_gap,1); % mean of the ay in the first half of the interval
%                 a_mean_secondhalf = zeros(last_index-min_gap,1); % mean of the ay in the second half of the interval
% %                 rdot_mean = zeros(last_index-min_gap,1); % Mean of r within each interval
% %                 rdot_mean_firsthalf = zeros(last_index-min_gap,1); % Mean of r within each interval
% %                 rdot_mean_secondhalf = zeros(last_index-min_gap,1); % Mean of r within each interval
%                 for ct1=min_gap+1:last_index
%                     coeff(ct1-min_gap,:) = [ones(ct1,1) y(1:ct1)]\r(1:ct1);
%                     r_mean(ct1-min_gap) = mean(r(1:ct1));
%                     a_mean(ct1-min_gap) = mean(a(1:ct1));
% %                     rdot_mean(ct1-min_gap) = mean(rdot(1:ct1));
%                     if rem(ct1,2) == 0
% %                         rdot_mean_firsthalf(ct1-min_gap) = mean(rdot(1:ct1/2));
% %                         rdot_mean_secondhalf(ct1-min_gap) = mean(rdot(ct1/2+1:ct1));
% 
%                         a_mean_firsthalf(ct1-min_gap) = mean(a(1:ct1/2));
%                         a_mean_secondhalf(ct1-min_gap) = mean(a(ct1/2+1:ct1));
%                         
%                         coeff_firsthalf(ct1-min_gap,:) = [ones(ct1/2,1) y(1:ct1/2)]\r(1:ct1/2);
%                         coeff_secondhalf(ct1-min_gap,:) = [ones(ct1/2,1) y(ct1/2+1:ct1)]\r(ct1/2+1:ct1);
%                     else
% %                         rdot_mean_firsthalf(ct1-min_gap) = mean(rdot(1:floor(ct1/2)+1));
% %                         rdot_mean_secondhalf(ct1-min_gap) = mean(rdot(floor(ct1/2)+1:ct1));
%                         a_mean_firsthalf(ct1-min_gap) = mean(a(1:floor(ct1/2)+1));
%                         a_mean_secondhalf(ct1-min_gap) = mean(a(floor(ct1/2)+1:ct1));  
% 
%                         coeff_firsthalf(ct1-min_gap,:) = [ones(floor(ct1/2)+1,1) y(1:floor(ct1/2)+1)]\r(1:floor(ct1/2)+1);
%                         coeff_secondhalf(ct1-min_gap,:) = [ones(floor(ct1/2)+1,1) y(floor(ct1/2)+1:ct1)]\r(floor(ct1/2)+1:ct1); 
%                     end
%                 end
%                 
%                 for ct_factor=1:length(factors)
%                     factor = factors(ct_factor);
%                     
%                     indx = find(abs(coeff(:,1)-r_mean) <= factor*params(1) & abs(coeff(:,2)) <= factor*params(2) & ...
%                                 abs(coeff_firsthalf(:,1)-r_mean) <= factor*2*params(1) & abs(coeff_firsthalf(:,2)) <= factor*2*params(2) & ...
%                                 abs(coeff_secondhalf(:,1)-r_mean) <= factor*2*params(1) & abs(coeff_secondhalf(:,2)) <= factor*2*params(2) & ...
%                                 a_mean <= 0 & a_mean >= -4 & ...
%                                 a_mean_firsthalf <= 0 & a_mean_firsthalf >= -4 & ...
%                                 a_mean_secondhalf <= 0 & a_mean_secondhalf >= -4);
%                     if ~isempty(indx)
%                         possible_indices{ct_factor} = [possible_indices{ct_factor}; ct*ones(numel(indx),1) ct+min_gap+indx-1];
%                     end
%                     
%                 end
%                 
%             end
%             
%             for ct_factor=1:length(factors)
%                 rmse_constFit = arrayfun(@(idx) rmse_func(idx, possible_indices{ct_factor}, r_all), 1:size(possible_indices{ct_factor},1));
%                 if ~isempty(possible_indices{ct_factor})
%                     intervalArray = rrefEstimateInterval.createIntervalArray(possible_indices{ct_factor}, rmse_constFit');
%                     
%                     % Time-window dependent best intervals for a particular
%                     % factor
%                     isChosenInterval = rrefEstimateInterval.findRepresentativeIntervalsPerTimeWindow(intervalArray);
%                     intervalArray = intervalArray(isChosenInterval);
%                     intervals_td = possible_indices{ct_factor}(isChosenInterval,:);
%                     
%                     [intervalArray.factor] = deal(factors{ct_factor});
%                     
%                     values = num2cell(arrayfun(@(i,j) sum(y_all(i:j).*v_all(i:j))/sum(y_all(i:j).^2),intervals_td(:,1),intervals_td(:,2)));
%                     [intervalArray.rref] = values{:};
%                     
%                     values = num2cell(arrayfun(@(i,j)  mean(y_all(i:j)),intervals_td(:,1),intervals_td(:,2)));
%                     [intervalArray.ymean] = values{:};
%                     
%                     values = num2cell(arrayfun(@(i,j)  mean(v_all(i:j)),intervals_td(:,1),intervals_td(:,2)));
%                     [intervalArray.vmean] = values{:};
%                     
%                     values = num2cell(arrayfun(@(i,j)  mean(r_all(i:j)),intervals_td(:,1),intervals_td(:,2)));
%                     [intervalArray.rmean] = values{:};
%                     
%                     obj.rrefSegments_tw_fac = [obj.rrefSegments_tw_fac; intervalArray];
%                     
%                     
%                     
%                     % Time-window independent best intervals
%                     chosenIntervals_ti = Interval.findRepresentativeIntervals(intervalArray);
%                     intervals_ti = possible_indices{ct_factor}(chosenIntervals_ti,:);
%                     obj.rrefSegments(ct_factor).intervals_ti = intervals_ti;
%                     obj.rrefSegments(ct_factor).rref_ti = arrayfun(@(i,j) sum(y_all(i:j).*v_all(i:j))/sum(y_all(i:j).^2),intervals_ti(:,1),intervals_ti(:,2));
%                     obj.rrefSegments(ct_factor).ymean_ti = arrayfun(@(i,j) mean(y_all(i:j)),intervals_ti(:,1),intervals_ti(:,2));
%                     obj.rrefSegments(ct_factor).vmean_ti = arrayfun(@(i,j) mean(v_all(i:j)),intervals_ti(:,1),intervals_ti(:,2));
%                     obj.rrefSegments(ct_factor).rmean_ti = arrayfun(@(i,j) mean(r_all(i:j)),intervals_ti(:,1),intervals_ti(:,2));
%                     coeff = arrayfun(@(i,j) [[ones(j-i+1,1) y_all(i:j)]\r_all(i:j)]',intervals_ti(:,1),intervals_ti(:,2),'UniformOutput',false);
%                     coeff = vertcat(coeff{:});
%                     obj.rrefSegments(ct_factor).const_rvsy_ti = coeff(:,1);
%                     obj.rrefSegments(ct_factor).slope_rvsy_ti = coeff(:,2);
%                     
%                     obj.rrefSegments(ct_factor).xTravelled_ti = arrayfun(@(i,j) abs(max(x_all(i:j)) - min(x_all(i:j))), intervals_ti(:,1), intervals_ti(:,2)) ;
%                     obj.rrefSegments(ct_factor).yTravelled_ti = arrayfun(@(i,j) abs(max(y_all(i:j)) - min(y_all(i:j))), intervals_ti(:,1), intervals_ti(:,2)) ;
%                     obj.rrefSegments(ct_factor).zTravelled_ti = arrayfun(@(i,j) abs(max(z_all(i:j)) - min(z_all(i:j))), intervals_ti(:,1), intervals_ti(:,2)) ;
%                     obj.rrefSegments(ct_factor).xmean_ti = arrayfun(@(i,j) mean(x_all(i:j)), intervals_ti(:,1), intervals_ti(:,2)) ;
%                     obj.rrefSegments(ct_factor).zmean_ti = arrayfun(@(i,j) mean(z_all(i:j)), intervals_ti(:,1), intervals_ti(:,2)) ;
%                     obj.rrefSegments(ct_factor).fd_analytical_ti = arrayfun(@(i,j) log(abs(y_all(j)/y_all(i))), intervals_ti(:,1), intervals_ti(:,2))./obj.rrefSegments(ct_factor).rref_ti; % analytically computed flight duration 
%                     obj.rrefSegments(ct_factor).fd_actual_ti = arrayfun(@(i,j) diff(t_all([i j])), intervals_ti(:,1), intervals_ti(:,2)); % actual flight duration
%                     
%                     % Time-window dependent best intervals for a particular
%                     % factor
%                     chosenIntervals_td = Interval.findRepresentativeIntervalsPerTimeWindow(intervalArray);
%                     intervals_td = possible_indices{ct_factor}(chosenIntervals_td,:);
%                     obj.rrefSegments(ct_factor).intervals_td = intervals_td;
%                     obj.rrefSegments(ct_factor).rref_td = arrayfun(@(i,j) sum(y_all(i:j).*v_all(i:j))/sum(y_all(i:j).^2),intervals_td(:,1),intervals_td(:,2));
%                     obj.rrefSegments(ct_factor).ymean_td = arrayfun(@(i,j) mean(y_all(i:j)),intervals_td(:,1),intervals_td(:,2));
%                     obj.rrefSegments(ct_factor).vmean_td = arrayfun(@(i,j) mean(v_all(i:j)),intervals_td(:,1),intervals_td(:,2));
% %                     
%                     
%                 end
%             end
%             
%             % Throw all intervals (different factors and time-windows) together and find best representatives
%             
%             
%         end
        
    end
    
    methods(Static)
        function in = IsInsideCylinder(p1, p2, radius, q)
            % p1 - 1X3
            % p2 - 1X3
            % r - 1X1
            % q - query points (NX3)
            vec = p2-p1;
            const = radius*norm(vec);
            
            in = dot(q-p1, repmat(vec,size(q,1),1), 2) >= 0 & ...
                dot(q-p2, repmat(vec,size(q,1),1), 2) <= 0 & ...
                sum((cross(q-p1, repmat(vec, size(q,1), 1))).^2,2).^0.5 <= const;
            
        end
        
        function dydt = const_drdt_movement(t, y, rdot)
            % Differential equations governing flying at constant dr/dt 
            % used in ode45 for integration
            % rdot - the constant value of dr/dt at which bbee is flying
            
            dydt(1,1) = -y(2);
            dydt(2,1) = (rdot*y(1)^2-y(2)^2)./y(1);
            
        end
    end
end