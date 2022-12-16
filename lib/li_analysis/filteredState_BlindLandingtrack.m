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
                    any(obj.filteredState(:,4) < -0.18)
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
            ylim([0 ceil(-1*ceil(min(data4rref(:,6)./data4rref(:,3))-1))+1]);
            yticks(0:1:ceil(-1*ceil(min(data4rref(:,6)./data4rref(:,3))-1))+1);
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
        
        function plotHandles = plot_rrefs_with_rawData(obj, factor, raw_state)
            % Plots V vs y and V/y vs y highling change in r* within the same track
            % excerpt with positive V and positive y
            % This function plots time-window independent best intervals
            
            % factor - for which threshold factor plots are required to be
            % produced
            
            % This function also overlays rawData
            % rawState has 18 columns (obj_id	frame	timestamp	x	y	z	xvel	yvel	zvel P00	P01	P02	P11	P12	P22	P33	P44	P55)
        
            
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
            complete_state = obj.filteredState(indx_start:end,:); % Butterworth filtered state - 10 columns (timestamp	x	y	z	xvel	yvel	zvel xacc yacc zacc)
            
            plotHandles = figure;
            p1 = subplot(2,1,1); hold on;
%             plot(-raw_state(:,5),raw_state(:,8),'.','MarkerSize',13,'MarkerFaceColor',[154,205,50]./255, 'MarkerEdgeColor',[154,205,50]./255');
            
            plot(-complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-raw_state(:,5),raw_state(:,8),'k');
            ylabel('V (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3)], 'FontSize', 16);
            %     ylim([-0.2 p1.YLim(2)]);
            %     ylim([-0.2 0.6]);
            %     yticks([-0.2:0.2:0.6]);
            %     xticks([p1.XLim(1):0.1:0]);
            
            
            p2 = subplot(2,1,2); hold on;
            %     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
%             plot(-raw_state(:,5),raw_state(:,8)./-raw_state(:,5),'.','MarkerSize',13,'MarkerFaceColor',[154,205,50]./255, 'MarkerEdgeColor',[154,205,50]./255');
            
            plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-raw_state(:,5),raw_state(:,8)./-raw_state(:,5),'k');
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
            elseif strcmpi(obj.landingSide, 'Feeder')
                [obj.rrefSegments.side] = deal(2);
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
        
%         function find_rrefEntry(obj)
%             % Goes through obj.rrefSegments and finds track segments
%             % corresponding to entry into constant r.
%             % monotonic increase or decrease (decided by first point just
%             % before start of constant-r) decides the start of entry
%             
%             if isempty(obj.rrefSegments)
%                 return
%             end
%             
%             indx = size(obj.filteredState,1);
%             t_all = obj.filteredState(1:indx,1)-obj.filteredState(1,1);
%             y_all = obj.filteredState(1:indx,3);
%             v_all = obj.filteredState(1:indx,6);
%             a_all = obj.filteredState(1:indx,9);
%             r = v_all./y_all;
%             diffr = diff(r);
%             
% %             figure;
% %             plot(y_all,r,'o');
% %             text(y_all,r,string(1:length(r)));
% %             title('y vs r')
%             
%             for ct=1:length(obj.rrefSegments)
%                 obj.rrefEntrySegments(ct) = data4rrefEntry();
%                 
%                 % Saving some parameters
%                 obj.rrefEntrySegments(ct).factor = obj.rrefSegments(ct).factor;
%                 
%                 % Find entry segment
%                 if isempty(obj.rrefSegments(ct).intervals_ti)
%                     continue
%                 end
%                 
%                 
%                 rref_intervals = [1 1; sortrows(obj.rrefSegments(ct).intervals_ti,1,'ascend')];
%                 indices = nan(size(rref_intervals,1)-1,1);
%                 for ct1=2:size(rref_intervals,1)
%                     if obj.rrefSegments(ct).rmean_ti(ct1-1) > -0.5
%                         continue
%                     end
%                     
%                     if diffr(rref_intervals(ct1,1)) > 0
%                         indx1 = find(diffr(rref_intervals(ct1-1,2):rref_intervals(ct1,1)) < 0, 1, 'last');
%                     elseif diffr(rref_intervals(ct1,1)) < 0
%                         indx1 = find(diffr(rref_intervals(ct1-1,2):rref_intervals(ct1,1)) > 0, 1, 'last');
%                     end
%                     if isempty(indx1)
%                         indx1 = rref_intervals(ct1-1,2);
%                     else
%                         indx1 = indx1 + rref_intervals(ct1-1,2);
%                     end
%                     
%                     indx2 = find(r(rref_intervals(ct1-1,2):rref_intervals(ct1,1)) > -0.5, 1, 'last');
%                     if isempty(indx2)
%                         indx2 = rref_intervals(ct1-1,2);
%                     else
%                         indx2 = indx2 + rref_intervals(ct1-1,2);
%                     end
%                     
%                     indx = max([indx1, indx2]);
%                     
%                     if rref_intervals(ct1,1)-indx>=15 && (obj.rrefSegments(ct).rmean_ti(ct1-1)*1.1 >= mean(r(indx:rref_intervals(ct1,1))) || ...
%                                                           obj.rrefSegments(ct).rmean_ti(ct1-1)*0.9 <= mean(r(indx:rref_intervals(ct1,1))))
%                         indices(ct1-1) = indx;
%                     end
%                 end
%                 
%                 if size(rref_intervals,1)>1
%                     [rref_intervals, sort_indices] = sortrows(obj.rrefSegments(ct).intervals_ti,1,'ascend');
%                     nonnan_sort_indices = sort_indices(~isnan(indices));
%                     obj.rrefEntrySegments(ct).intervals = [indices(~isnan(indices)) rref_intervals(~isnan(indices),:)]; % N by 2, where N is the number of non-overlapping time-window independent intervals
%                     obj.rrefEntrySegments(ct).rref = obj.rrefSegments(ct).rref_ti(nonnan_sort_indices);
%                     obj.rrefEntrySegments(ct).rmean = obj.rrefSegments(ct).rmean_ti(nonnan_sort_indices);
%                 end
%             end
%             
%             
%         end
%         
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
                    obj.rrefEntrySegments(ct).delta_Ventry(ct1,1) = diff(v_interval([1 end]));
                    obj.rrefEntrySegments(ct).delta_tentry(ct1,1) = diff(t_interval([1 end]));
                    obj.rrefEntrySegments(ct).amean_entry(ct1,1) = mean(a_interval);
                end
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
        
        function plotHandles = plot_rrefs_parallax(obj, factor)
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
            p1 = subplot(4,1,1); hold on;
            plot(-complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('V (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3)], 'FontSize', 16);
            %     ylim([-0.2 p1.YLim(2)]);
            %     ylim([-0.2 0.6]);
            %     yticks([-0.2:0.2:0.6]);
            %     xticks([p1.XLim(1):0.1:0]);
            
            
            p2 = subplot(4,1,2); hold on;
            %     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylabel('r (1/s)', 'FontSize', 16);
            %     ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
%             xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            ylim([0 ceil(-1*ceil(min(data4rref(:,6)./data4rref(:,3))-1))]);
            yticks(0:1:ceil(-1*ceil(min(data4rref(:,6)./data4rref(:,3))-1)));
            %     xticks(p1.XLim(1):0.1:0);
            
            p3 = subplot(4,1,3); hold on;
            plot(-complete_state(:,3),complete_state(:,5)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylim([-3 3]);
            ylabel('Horizontal flow', 'FontSize', 16);
            set(gca, 'FontSize', 16); %grid on;
            
            p4 = subplot(4,1,4); hold on;
            plot(-complete_state(:,3),complete_state(:,7)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            ylim([-3 3]);
            ylabel('Vertical flow', 'FontSize', 16);
            xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 16);
            
            linkaxes([p1, p2, p3, p4],'x');
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
                subplot(4,1,1);
                plot(-state_subset(:,3),state_subset(:,6),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                plot(-state_subset(:,3),-obj.rref*-state_subset(:,3),'LineWidth',2,'Color',[69 117 180]./255);
                plot([0 -state_subset(end,3)],[0 -obj.rref*-state_subset(end,3)],'--','LineWidth',2,'Color',[69 117 180]./255);
                
                subplot(4,1,2);
                plot(-state_subset(:,3),state_subset(:,6)./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
                %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
                plot(-state_subset([1,end],3),[-obj.rref, -obj.rref],'LineWidth',2,'Color',[69 117 180]./255);
                plot([p1.XLim(1) -state_subset(1,3)],[-obj.rref -obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
            end
            
        end
        
        function [plotHandles, AIC] = compute_and_compare_three_landing_strategies(obj, factor, createPlots)
            % This function computes and compares tracks for 3 landing
            % stratgies: const_r, const_taudot, step_rref
            % This is only done if there are more than one rref segment
            % identified for factor (an input) f.
            % For const_r, it is assumed that bee continues to fly the
            % first identified rref
            % For const_taudot, tractory is computed based on taudot (an
            % input to the function). It is usually the average taudot
            % observed in whole dataset.
            
            % Comparison is done based on AIC and is done only for the
            % parts of the track where rref segments are identified.
            
            
            % Plots V vs y and r vs y highlighting three strategies
            
            % factor - for which threshold factor plots are required to be
            % produced
            
            assert(~isempty(obj.rrefSegments) && length(factor) == 1);
            ct_factor = find(abs([obj.rrefSegments.factor] - factor) < 1e-6);
            if isempty(ct_factor)
                error('Can NOT find the r* intervals for the asked factor.');
            end
            intervals = obj.rrefSegments(ct_factor).intervals_ti;
            if isempty(intervals) || size(intervals,1) < 2 % only for tracks containing more than one rref segment
                plotHandles = [];
                dataout = [];
                return;
            end
            
            dataPerTrackExcerpt = struct.empty;
            for ct=1:size(intervals,1)
                dataPerTrackExcerpt(ct).actualstate = [-obj.filteredState(intervals(ct,1):intervals(ct,2),3) obj.filteredState(intervals(ct,1):intervals(ct,2),6) obj.filteredState(intervals(ct,1):intervals(ct,2),1)]; % [y V t]
                dataPerTrackExcerpt(ct).rref = obj.rrefSegments(ct_factor).rref_ti(ct);
%                 dataPerTrackExcerpt(ct).state4const_r = [t(intervals(ct,1)-intervals(ct,1)+1:intervals(ct,2)-intervals(ct,1)+1,:) state_const_r(intervals(ct,1)-intervals(ct,1)+1:intervals(ct,2)-intervals(ct,1)+1,:)];
%                 dataPerTrackExcerpt(ct).state4const_taudot = [t(intervals(ct,1)-intervals(ct,1)+1:intervals(ct,2)-intervals(ct,1)+1,:) state_const_taudot(intervals(ct,1)-intervals(ct,1)+1:intervals(ct,2)-intervals(ct,1)+1,:)];
%                 dataPerTrackExcerpt(ct).state_actual = obj.filteredState(intervals(ct,1):intervals(ct,2),[1 3 6]);
                dataPerTrackExcerpt(ct).state4step_rref = [-obj.filteredState(intervals(ct,1):intervals(ct,2),3) obj.rrefSegments(ct_factor).rref_ti(ct)*obj.filteredState(intervals(ct,1):intervals(ct,2),3) obj.filteredState(intervals(ct,1):intervals(ct,2),1)];
            end
            dummy_state = vertcat(dataPerTrackExcerpt.state4step_rref);
            [~,indx,~] = unique(dummy_state(:,3)); % column 3 is time
            state4step_rref = dummy_state(indx,:);
            
            dummy_state = vertcat(dataPerTrackExcerpt.actualstate);
            [~,indx,~] = unique(dummy_state(:,3)); % column 3 is time
            actualstate = dummy_state(indx,:);
            
            
            rrefSegs = obj.rrefSegments(ct_factor);
%             time_vec = obj.filteredState(rrefSegs.intervals_ti(1,1):end,1);
%             y_vec = -obj.filteredState(rrefSegs.intervals_ti(1,1):end,3);
%             v_vec = obj.filteredState(rrefSegs.intervals_ti(1,1):end,6);
%             a_vec = obj.filteredState(rrefSegs.intervals_ti(1,1):end,9);
%             
            time_vec = obj.filteredState(rrefSegs.intervals_ti(1,1):rrefSegs.intervals_ti(end,2),1);
            y_vec = -obj.filteredState(rrefSegs.intervals_ti(1,1):rrefSegs.intervals_ti(end,2),3);
            v_vec = obj.filteredState(rrefSegs.intervals_ti(1,1):rrefSegs.intervals_ti(end,2),6);
            a_vec = obj.filteredState(rrefSegs.intervals_ti(1,1):rrefSegs.intervals_ti(end,2),9);
            r_vec = v_vec./y_vec;
            taudot_vec = -1-y_vec.*a_vec./v_vec.^2;
            
            
            % simulate const_rref movement
            rref = mean(r_vec); % -rrefSegs.rref_ti(1); % assuming bee flies at first identified rref
            initialCondition = [-obj.filteredState(rrefSegs.intervals_ti(1,1),3) -rref*obj.filteredState(rrefSegs.intervals_ti(1,1),3)];
            [t_const_r,state_const_r] = ode45(@(t,y) filteredState_BlindLandingtrack.const_r_movement(t,y,rref), time_vec, initialCondition);

            % simulate const_taudot movement
            % assuming bee flies at the average taudot/m identified for whole dataset
            taudot = mean(taudot_vec);
            initialCondition = [-obj.filteredState(rrefSegs.intervals_ti(1,1),3) obj.filteredState(rrefSegs.intervals_ti(1,1),6)];
            [t_const_taudot,state_const_taudot] = ode45(@(t,y) filteredState_BlindLandingtrack.const_taudot_movement(t,y,taudot), time_vec, initialCondition);
            
            in_rref = false(length(time_vec),1);
            for ct=1:size(intervals,1)
                in_rref(rrefSegs.intervals_ti(ct,1)-rrefSegs.intervals_ti(1,1)+1:rrefSegs.intervals_ti(ct,2)-rrefSegs.intervals_ti(1,1)+1) = true;
            end
%             sum([state4step_rref(:,3)-t_const_r(in_rref)]) % for testing
            [actualstate(1,3)-t_const_r(1)]
            output = 2; % 1:2 for both y and V as output
            AIC.step_rref = filteredState_BlindLandingtrack.computeAIC([actualstate(:,output)-state4step_rref(:,output)]', size(intervals,1));
            AIC.const_r = filteredState_BlindLandingtrack.computeAIC([actualstate(:,output)-state_const_r(in_rref,output)]', 1);
            AIC.const_taudot = filteredState_BlindLandingtrack.computeAIC([actualstate(:,output)-state_const_taudot(in_rref,output)]', 1);
            
            if ~createPlots
                plotHandles = [];
                return;
            end
            
            tend = obj.filteredState(end,1);
            t_all = obj.filteredState(1:end,1);
            y_all = -obj.filteredState(1:end,3);
            v_all = obj.filteredState(1:end,6);
            a_all = obj.filteredState(1:end,9);
            r_all = v_all./y_all;
            taudot_all = -1-y_all.*a_all./v_all.^2;
            
            plotHandles = figure;
            figure(plotHandles);
                
            subplotHandles(1) = subplot(4,1,1); hold on;
            plot(t_all-tend, y_all,'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(t_const_r-tend, state_const_r(:,1),'LineWidth',3,'Color',[69 117 180]./255);
            plot(t_const_taudot-tend, state_const_taudot(:,1),'LineWidth',2,'Color',[255 0 0]./255);
            rrefState = {};
            for ct=1:size(intervals,1)
                rrefState{ct} = obj.filteredState(intervals(ct,1):intervals(ct,2),:);
                plot(rrefState{ct}(:,1)-tend, -rrefState{ct}(:,3),'LineWidth',2,'Color',[0 0 0]./255); 
                patch([rrefState{ct}(1,1) rrefState{ct}(1,1) rrefState{ct}(end,1) rrefState{ct}(end,1)]-tend, ...
                    [subplotHandles(1).YLim(1)  subplotHandles(1).YLim(2) subplotHandles(1).YLim(2)  subplotHandles(1).YLim(1)], ...
                    [.7 .7 .7], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            end
            ylabel('y (m)', 'FontSize', 15);
%             legend({'Actual data', 'constant-r approach', 'constant-\dot{\tau} approach', 'hybrid'}, 'FontSize', 13);
            set(gca, 'FontSize', 15);

            subplotHandles(2) = subplot(4,1,2); hold on;
            plot(t_all-tend, v_all,'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(t_const_r-tend, state_const_r(:,2),'LineWidth',3,'Color',[69 117 180]./255);
            plot(t_const_taudot-tend, state_const_taudot(:,2),'LineWidth',2,'Color',[255 0 0]./255);
            rrefState = {};
            for ct=1:size(intervals,1)
                rrefState{ct} = obj.filteredState(intervals(ct,1):intervals(ct,2),:);
                rref = obj.rrefSegments.rref_ti(ct);
                patch([rrefState{ct}(1,1) rrefState{ct}(1,1) rrefState{ct}(end,1) rrefState{ct}(end,1)]-tend, ...
                    [subplotHandles(2).YLim(1)  subplotHandles(2).YLim(2) subplotHandles(2).YLim(2)  subplotHandles(2).YLim(1)], ...
                    [.7 .7 .7], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
                plot(rrefState{ct}(:,1)-tend, rref.*rrefState{ct}(:,3),'LineWidth',2,'Color',[0 0 0]./255);
            end
            ylabel('V (m/s)', 'FontSize', 15);
            set(gca, 'FontSize', 15);

            subplotHandles(3) = subplot(4,1,3); hold on;
            plot(t_all-tend, r_all,'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(t_const_r-tend, state_const_r(:,2)./state_const_r(:,1),'LineWidth',3,'Color',[69 117 180]./255);
            plot(t_const_taudot-tend, state_const_taudot(:,2)./state_const_taudot(:,1),'LineWidth',2,'Color',[255 0 0]./255);
            rrefState = {};
            for ct=1:size(intervals,1)
                rrefState{ct} = obj.filteredState(intervals(ct,1):intervals(ct,2),:);
                rref = obj.rrefSegments.rref_ti(ct);
                patch([rrefState{ct}(1,1) rrefState{ct}(1,1) rrefState{ct}(end,1) rrefState{ct}(end,1)]-tend, ...
                    [subplotHandles(3).YLim(1)  subplotHandles(3).YLim(2) subplotHandles(3).YLim(2)  subplotHandles(3).YLim(1)], ...
                    [.7 .7 .7], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
                plot(rrefState{ct}(:,1)-tend, -rref*ones(length(rrefState{ct}(:,1)),1),'LineWidth',2,'Color',[0 0 0]./255);
            end
            ylabel('r (1/s)', 'FontSize', 15);
            set(gca, 'FontSize', 15);
            
            subplotHandles(4) = subplot(4,1,4); hold on;
            plot(t_all-tend, taudot_all,'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
            plot(t_const_r-tend, zeros(length(t_const_r),1),'LineWidth',6,'Color',[69 117 180]./255);
            plot(t_const_taudot-tend, taudot*ones(length(t_const_taudot),1),'LineWidth',2,'Color',[255 0 0]./255);
            rrefState = {};
            for ct=1:size(intervals,1)
                rrefState{ct} = obj.filteredState(intervals(ct,1):intervals(ct,2),:);
                patch([rrefState{ct}(1,1) rrefState{ct}(1,1) rrefState{ct}(end,1) rrefState{ct}(end,1)]-tend, ...
                    [subplotHandles(4).YLim(1)  subplotHandles(4).YLim(2) subplotHandles(4).YLim(2)  subplotHandles(4).YLim(1)], ...
                    [.7 .7 .7], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
                plot(rrefState{ct}(:,1)-tend, zeros(size(rrefState{ct},1),1),'LineWidth',2,'Color',[0 0 0]./255);
            end
            ylabel('\dot{\tau}', 'FontSize', 15,'interpreter','latex');
            xlabel('t (s)', 'FontSize', 15);
            set(gca, 'FontSize', 15);
                
            linkaxes(subplotHandles,'x');
            
            rrefState = vertcat(rrefState{:});
            tlim = [floor(min(rrefState(:,1)-tend)*10) ceil(max(rrefState(:,1)-tend)*10)]./10;
            for ct=1:4
                subplotHandles(ct).XLim = tlim;
                if ct==1
                    ylim = [0 ceil(max(-rrefState(:,3))*10)./10+0.1];
                elseif ct==2
                    ylim = [0 ceil(max(rrefState(:,6))*10)./10+0.2];
                elseif ct==3
                    r = -rrefState(:,6)./rrefState(:,3);
                    ylim = [floor(min(r)*10)./10-1 ceil(max(r)*10)./10+1];
                elseif ct==4
                    dummy = taudot_all(t_all-tend>=tlim(1) & t_all-tend<=tlim(2));
                    ylim = [floor(min(dummy)*10)./10-0.2 ceil(max(r)*10)./10+0.2];
                end
                subplotHandles(ct).YLim = ylim;
            end
            
            
    end
        
        function [Vhybrid, Vconst_r] = compute_avgV(obj, factor)
            assert(~isempty(obj.rrefSegments) && length(factor) == 1);
            ct_factor = find(abs([obj.rrefSegments.factor] - factor) < 1e-6);
            if isempty(ct_factor)
                error('Can NOT find the r* intervals for the asked factor.');
            end
            intervals = obj.rrefSegments(ct_factor).intervals_ti;
            if size(intervals,1) < 2 % only for tracks containing more than one rref segment
                Vhybrid = [];
                Vconst_r = [];
                return;
            end
            
            Vhybrid = (obj.filteredState(intervals(1,1),3)-obj.filteredState(intervals(end,2),3))./...
                (obj.filteredState(intervals(1,1),1)-obj.filteredState(intervals(end,2),1));
%             Vhybrid = mean(obj.filteredState(intervals(1,1):intervals(end,2),6));
            Vconst_r = obj.rrefSegments(ct_factor).rref_ti(1)*mean(obj.filteredState(intervals(1,1):intervals(end,2),3));
%             ystart = -obj.filteredState(intervals(1,1),3);
%             yend = -obj.filteredState(intervals(end,2),3);
%             Vconst_r = 2/(-obj.rrefSegments(ct_factor).rref_ti(1)*(ystart+yend));
        end
        
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
            % rdot - the constant value of dr/dt at which bee is flying
            
            dydt(1,1) = -y(2);
            dydt(2,1) = (rdot*y(1)^2-y(2)^2)./y(1);
            
        end
        
        function dydt = const_r_movement(t,y,rref)
            % Differential equations governing flying at constant dr/dt 
            % used in ode45 for integration
            % rref - the constant value of r at which bee is flying
            dydt(1,1) = -y(2);
            dydt(2,1) = -rref*y(2);
        end
        
        function dydt = const_taudot_movement(t,y,taudot)
            % Differential equations governing flying at constant taudot 
            % used in ode45 for integration
            % taudot - the constant value of dr/dt at which bee is flying
            dydt(1,1) = -y(2);
            dydt(2,1) = -(taudot*y(2)^2+y(2)^2)./y(1);
            
        end
        
        function AIC = computeAIC(error, nparams)
            % error: ny by N matrix
            % ny: no. of outputs
            % N: number of values in the estimation dataset
            % nparams: number of estimated parameters
            
            N = size(error,2);
            ny = size(error,1);
            errorSum = 0;
            for ct=1:N
                errorSum = errorSum + error(:,ct)*error(:,ct)';
            end
            AIC = N*log(det(errorSum./N)) + 2*nparams + N*(ny*(log(2*pi)+1));
        end
        
    end
end