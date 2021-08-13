classdef filteredState_BlindLandingtrack < handle
    properties (Access = public)
        filteredState = [];  % Butterworth filtered state - 10 columns (timestamp	x	y	z	xvel	yvel	zvel xacc yacc zacc)
        landingSide = '';
        
        % to store parameters for r* segments (these r* segments are
        % automatically computed
        rrefSegments = data4rrefEstimate.empty;
        
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
        
        function is = isGoldenFlight(obj, landingDiscs)
            % Determines if the landing started from opposite disc
            % (referred to as golden flight)
            
            % landingDiscs - 1x2 LandingDisc array
            [~, ymin_indx] = min(obj.filteredState(:,3));
            state0 = obj.filteredState(ymin_indx,:);
            
            oppositeDisc = landingDiscs(arrayfun(@(x) ~strcmpi(x.side, obj.landingSide), landingDiscs));
            thisDisc = landingDiscs(arrayfun(@(x) strcmpi(x.side, obj.landingSide), landingDiscs));
            
            center = (oppositeDisc.center - thisDisc.center)';
            height_cylinder = 0.06;
            if filteredState_BlindLandingtrack.IsInsideCylinder(center - [0 height_cylinder/2 0], center + [0 height_cylinder/2 0], oppositeDisc.radius, state0(:,2:4))
                is = true;
            else
                is = false;
            end
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
        
        function compute_instabilityFollows(obj, yrange)
            % For every rref segment, checks whether instability follows
            % within yrange after it ends.
            % Condition for instability (V<0.05 m/s)
            
            vthreshold = 0.05; % in m/s
            
            if isempty(obj.rrefSegments)
                return
            end
            
            indx = size(obj.filteredState,1);
%             t_all = obj.filteredState(1:indx,1)-obj.filteredState(1,1);
            y_all = -obj.filteredState(1:indx,3);
            v_all = obj.filteredState(1:indx,6);
%             a_all = obj.filteredState(1:indx,9);
%             r = v_all./y_all;
            
            for ct=1:length(obj.rrefSegments)
                if isempty(obj.rrefSegments(ct).intervals_ti)
                    continue
                end
                
                rref_intervals = obj.rrefSegments(ct).intervals_ti;
                
                instabilityFollows = false(size(rref_intervals,1),1);
                for ct1=1:size(rref_intervals,1)
                    
                    
                    y_part = y_all(rref_intervals(ct1,2):end)-y_all(rref_intervals(ct1,2));
                    indx = find(y_part<-yrange,1);
                                        
                    if isempty(indx) 
                        continue;
                    elseif any(v_all(rref_intervals(ct1,2):rref_intervals(ct1,2)+indx-1)<vthreshold)
                        instabilityFollows(ct1) = true;
                    end
                    
                end
                obj.rrefSegments(ct).instabilityFollows = instabilityFollows;
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
    end
    
end