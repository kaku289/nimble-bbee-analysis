classdef BlindLandingtrack < handle
    properties (Access = public)
        
        foldername = ''; % file that contains the track
        
        datenum; % date on which the track was obtained
        pattern;
        patternnum;
        setID;
        beeID;
        flightID;
        filename;
        basedir;
        
        
        
        light = '';
        
        obj_id = []; 

        rawTrack = rawState_BlindLandingtrack.empty; % raw state, 16 columns, structure is created because # of landing tracks can vary for each object id
        state = filteredState_BlindLandingtrack.empty; % Butterworth filtered state - 10 columns (timestamp	x	y	z	xvel	yvel	zvel xacc yacc zacc)
       
        
        dayTimeInstant = ''; % number of seconds since start of the day
        treatmentTimeInstant = ''; % number of seconds since start of the treatment
                
        % parameters
        time; % time vector for plots 
        displacement; % Eucledian distance from 3D position of disc center
        
        vbyz = nan;
        
        dt = 0.00571966171264648; % delta t between two consecutive frames (in s) from hdf5 files
        
        % parameters in Landing Disc Frame (LDF)
        state_LDF = filteredState_BlindLandingtrack.empty; % N_excerpts by 1 vector
        raw_LDF = rawState_BlindLandingtrack.empty; % N_excerpts by 1 vector
               
        % To store parameters extracted using GUIDE
        DataGUI = DataGUI_BlindTracks.empty; % N_excerpts by 1 vector
  
        equinox_state_LDF = []; % point where first oscillation is observed in landing trajectory
        stable_state_LDF = []; % traj between equinox_state and some y e.g. 0.15m - To get info about r_ref
        
        % For estimating the parameters of the landing model
        track_subset_sysID = trackForLandingModel.empty; % N_excerpts by 1 vector
                
    end
    methods
        function obj = BlindLandingtrack(pattern, light, foldername, datenum, obj_id)
            if nargin == 5
                obj.pattern = pattern;
                obj.light = light;
                obj.foldername = foldername;
                obj.datenum = datenum;
                obj.obj_id = obj_id;
            elseif nargin == 0
            end
%             
%             if ~isempty(rawTrack)
%                 
%                 obj.confirmedChecks_le = confirmedChecks;
% 
%                 obj.rawTrack = rawTrack;
% 
%                 obj.obj_id = unique(rawTrack(:,1));
%                 obj.nObjects = length(obj.obj_id);
% 
%                 obj.scoredLandingParams = copy(videoParams.landingParams);
% 
%                 end_t = videoParams.timestamps(nanmax([obj.scoredLandingParams.legextension_framenumber ...
%                         obj.scoredLandingParams.touchdown_framenumber ...
%                         obj.scoredLandingParams.touchdown2_framenumber])) + 2*obj.dt;
% 
%                 obj.rawTrack = sortrows(obj.rawTrack(obj.rawTrack(:,3)-(end_t) <= 1e-5,:), [3 2]);
% 
%                 if ~isempty(obj.rawTrack) && sum(isnan(obj.rawTrack(:,3))) == 0
%                     obj.filterTrack(20);
%                 end
%             end
            
        end
        
        function state = filterTrack(obj, rawTrack, num, den)
            % rawTrack - instance of rawState_BlindLandingtrack
            % num - numerator for the butterworth filter tf
            % den - denominator for the butterworth filter tf
            
            % stateOutput - instance of filteredState_BlindLandingtrack
            
            
            % fs = 1/nanmean(diff(obj.rawTrack(:,3)));

            equidistant_time = [rawTrack.rawState(1,3):obj.dt:rawTrack.rawState(end,3)]';
            
            pos = interp1(rawTrack.rawState(:,3)-rawTrack.rawState(1,3), rawTrack.rawState(:,4:6), ...
                equidistant_time-equidistant_time(1), 'makima');
            
            filteredState = zeros(size(pos,1),10);
            filteredState(:,1) = equidistant_time;
            filteredState(:,2:4) = filtfilt(num, den, pos);
            filteredState(:,5:7) =  diffxy(equidistant_time, filteredState(:,2:4)); %filtfilt(num, den, diffxy(equidistant_time, obj.state(:,2:4)));
            filteredState(:,8:10) = diffxy(equidistant_time, filteredState(:,2:4), [], 2); %filtfilt(num, den, diffxy(equidistant_time, obj.state(:,2:4), [], 2));
            state = filteredState_BlindLandingtrack(filteredState, rawTrack.landingSide);
        end
        
        function filterRawTracks(obj, fc)
            fs = 1/obj.dt; % fs is 174.8355 Hz
            [num, den] = butter(2, fc/(fs/2),'low'); % 2nd order Butterworth filter %XXXX CHECK ORDER OF THE FILTER
      
            obj.state = filteredState_BlindLandingtrack.empty;
            for ct=1:length(obj.rawTrack)
                obj.state(end+1) = obj.filterTrack(obj.rawTrack(ct), num, den);
            end
            
        end
        function filterRawTracks_2(obj, fc) % Not being used
            fs = 1/obj.dt; % fs is 174.8355 Hz
            [num, den] = butter(2, fc/(fs/2),'low'); % 2nd order Butterworth filter %XXXX CHECK ORDER OF THE FILTER
      
            obj.state = filteredState_BlindLandingtrack.empty;
            toDelete = false(length(obj.rawTrack),1);
            for ct=1:length(obj.rawTrack)
                % For the tracks that has more than 10 time points gap, store
                % only the part that covers maximum delta y in between and
                % has at least 20 timepoints in the interval
                % 
                currentTrack = obj.rawTrack(ct);
                indices = find(diff(currentTrack.rawState(:,3))>10*obj.dt); % because we can't filter if the gap is more than 10 time points
                
                if ~isempty(indices)
                    indices1 = [0 indices];
                    indices2 = [indices1(2:end) size(currentTrack.rawState,1)];

                    yTravelled = arrayfun(@(i,j) abs(max(currentTrack.rawState(i:j,5))-min(currentTrack.rawState(i:j,5))), indices1+1, indices2);
                    nPoints = indices2-indices1;
                    dummy = sortrows([yTravelled' nPoints' indices1' indices2'], [-1 -2]);
                    indx = find(dummy(:,2)>20,1); % picking first interval in dummy that has at least 20 points in between them

                    if ~isempty(indx)
                        % % % save the state after filtering
                        % create rawState_BlindLandingtrack from the selected interval 
                        instance = currentTrack.createNewInstanceFromSubset(dummy(indx,3)+1,dummy(indx,4));
                        obj.state(end+1) = obj.filterTrack(instance, num, den);
                    else
                        % track becomes smaller than 20 time points after
                        % removal of time gaps, therefore delete it.
                        toDelete(ct) = true;
                    end
                else
                    % % % save the state after filtering
                    obj.state(end+1) = obj.filterTrack(obj.rawTrack(ct), num, den);
                end
                
            end
            obj.rawTrack(toDelete) = [];
        end
        
        function removeDuplicateTracks(obj)
            % This function removes the duplicates of a raw track that
            % might come from the same object oscillating around point of
            % reference (e.g., y = 0.10m plane)
            
            isDuplicate = false(length(obj.rawTrack),1);
            for ct=1:length(obj.rawTrack)-1
                currentTrack = obj.rawTrack(ct);
                for ct1=ct+1:length(obj.rawTrack)
                    if ~isDuplicate(ct1)
                        duplicate = intersect(currentTrack.rawState, obj.rawTrack(ct1).rawState, 'rows');
                        if size(duplicate,1)/size(currentTrack.rawState,1) > 0.8
                            isDuplicate(ct) = true;
                            break;
                        elseif size(duplicate,1)/size(obj.rawTrack(ct1).rawState,1) > 0.8
                            isDuplicate(ct1) = true; 
                        elseif size(duplicate,1)/size(currentTrack.rawState,1) > 0.4 && ...
                           size(duplicate,1)/size(obj.rawTrack(ct1).rawState,1) > 0.4 
                            isDuplicate(ct1) = true;                            
                        end
                    end
                end
            end
            obj.rawTrack(isDuplicate) = [];
            
        end
        
        function discardTrackPartsAfterGaps(obj)
            % This function removes the parts of the tracks that lie after
            % a gap of more than 10 time points
            % NOTE - Instead of deleting, the whole track (or its part) is kept in the
            % rawTrack and only the part that is away from the
            % platform is filtered and stored in state.
            toDelete = false(length(obj.rawTrack),1);
            for ct=1:length(obj.rawTrack)
                currentTrack = obj.rawTrack(ct);
                indx = find(diff(currentTrack.rawState(:,3))>10*obj.dt,1,'last');
                if ~isempty(indx) && size(currentTrack.rawState,1)-indx>20
                    obj.rawTrack(ct).rawState(1:indx,:) = [];
                elseif ~isempty(indx) && size(currentTrack.rawState,1)-indx<20
                    % track becomes smaller than 20 time points after
                    % removal of time gaps, therefore delete it.
                    toDelete(ct) = true;
                end
            end
            obj.rawTrack(toDelete) = [];
        end
        
        function compute_states_in_LDF(obj, landingDiscs)
            % Computes state in a inertial reference frame attached to the
            % center of the landing disc
            
            % landingDisc - array of instances of LandingDisc class
            if length(obj.rawTrack)~=length(obj.state) % if filtered data doesn't exist, filter it first
                obj.filterRawTracks(20);
            end
            
            obj.state_LDF = filteredState_BlindLandingtrack.empty;
            for ct=1:length(obj.state)
                landing_side = obj.state(ct).landingSide;
                landing_disc = landingDiscs(strcmpi({landingDiscs.side}, landing_side));
                
%                 assert(length(landing_disc) == 1, 'Object is landing at multiple discs! Not possible...');
                
                stateLDF = BlindLandingtrack.convert_to_landing_disc_reference_frame2(obj.state(ct).filteredState, landing_disc.center', landing_side);
                obj.state_LDF(end+1) = filteredState_BlindLandingtrack(stateLDF, landing_side);
            end
            
        end
        
        function compute_rawData_in_LDF(obj, landingDiscs)
            % Computes rawData in a inertial reference frame attached to the
            % center of the landing disc
            
            % landingDisc - array of instances of LandingDisc class
            
            obj.raw_LDF = rawState_BlindLandingtrack.empty;
            for ct=1:length(obj.rawTrack)
                landing_side = obj.rawTrack(ct).landingSide;
                landing_disc = landingDiscs(strcmpi({landingDiscs.side}, landing_side));
                
%                 assert(length(landing_disc) == 1, 'Object is landing at multiple discs! Not possible...');
                
                stateLDF = BlindLandingtrack.convert_to_landing_disc_reference_frame2(obj.rawTrack(ct).rawState(:,3:9), landing_disc.center', landing_side);
                obj.raw_LDF(end+1) = rawState_BlindLandingtrack([obj.rawTrack(ct).rawState(:,1:2) stateLDF obj.rawTrack(ct).rawState(:,10:end)], landing_side);
            end
            
        end
        
        function plotHandles = plotDataLDF_Time_rawVSfiltered(obj, ct_excerpt, subplotHandles)
            % Plot optical flow parameters (V, y, r with time)
            % These parameters are plotted with time
            
            % ct_excerpt - ct whose raw and filtered plots are required
            % subplotHandles - array of axes handles to plot data in
           

            filtered = obj.state_LDF(ct_excerpt).filteredState;
            raw = obj.raw_LDF(ct_excerpt).rawState;
            t0 = filtered(end,1);
            
            if isempty(subplotHandles)
                opticalExpansionPlot = figure;
                figure(opticalExpansionPlot);
                
                subplotHandles(1) = subplot(3,1,1); hold on;
                plot(filtered(:,1)-t0, filtered(:,6),'b.','MarkerSize',10);
                plot(raw(:,3)-t0, raw(:,8),'r.','MarkerSize',10);
                ylabel('V_{gy} (m/s)', 'FontSize', 15);
                legend({'filtered','raw'},'Location','best');
                set(gca, 'FontSize', 15); grid on;
%                 ylim([-8 2]);

                subplotHandles(2) = subplot(3,1,2); hold on;
                plot(filtered(:,1)-t0,...
                    filtered(:,3),'b.','MarkerSize',10);
                plot(raw(:,3)-t0, raw(:,5),'r.','MarkerSize',10);
                ylabel('y (m)', 'FontSize', 15);
    %             xlabel('Time (s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;

                subplotHandles(3) = subplot(3,1,3); hold on;
                plot(filtered(:,1)-t0,...
                    filtered(:,6)./filtered(:,3),'b.','MarkerSize',10);
                plot(raw(:,3)-t0, raw(:,8)./raw(:,5),'r.','MarkerSize',10);
                ylabel('r (1/s)', 'FontSize', 15);
                xlabel('Time (s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                
                plotHandles = opticalExpansionPlot;
                
            else
                
                assert(length(subplotHandles) == 3, 'Three axes handles are needed to plot the data');
                
%                 set(0,'CurrentFigure',subplotHandles(1).Parent);
%                 axes(subplotHandles(1));
                plot(filtered(:,1)-t0, filtered(:,6),'b.','MarkerSize',10, 'Parent', subplotHandles(1));
                plot(raw(:,3)-t0, raw(:,8),'r.','MarkerSize',10, 'Parent', subplotHandles(1));
                legend(subplotHandles(1),{'filtered','raw'},'Location','best');
%                 ylim([-8 2]);
                
%                 axes(subplotHandles(2));
                plot(filtered(:,1)-t0, filtered(:,3),'b.','MarkerSize',10, 'Parent', subplotHandles(2));
                plot(raw(:,3)-t0, raw(:,5),'r.','MarkerSize',10, 'Parent', subplotHandles(2));
%                 ylim([-60 60]);
                
%                 axes(subplotHandles(3));
                plot(filtered(:,1)-t0, filtered(:,6)./filtered(:,3),'b.','MarkerSize',10, 'Parent', subplotHandles(3));
                plot(raw(:,3)-t0, raw(:,8)./raw(:,5),'r.','MarkerSize',10, 'Parent', subplotHandles(3));
                
                plotHandles = [];
            end
        end
        
        
        function plotHandles = plotDataLDF_Distance_rawVSfiltered(obj, ct_excerpt, subplotHandles)
            % Plot optical flow parameters (V and r with y)
            % These parameters are plotted with distance from the platform
            
            % ct_excerpt - ct whose raw and filtered plots are required
            % subplotHandles - array of axes handles to plot data in
            
%             str = '#375E98';
%             color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
            assert(length(obj.state) == length(obj.rawTrack));
            
            filtered = obj.state_LDF(ct_excerpt).filteredState;
            raw = obj.raw_LDF(ct_excerpt).rawState;
            
            
            if isempty(subplotHandles)
                opticalExpansionPlot = figure;
                figure(opticalExpansionPlot);
                subplot(2,1,1); hold on;
                plot(filtered(:,3), filtered(:,6),'b.','MarkerSize',10);
                plot(raw(:,5), raw(:,8),'r.','MarkerSize',10);
                ylabel('V_{gy} (m/s)', 'FontSize', 15);
                legend({'filtered','raw'},'Location','best');
                set(gca, 'FontSize', 15); grid on;

                subplot(2,1,2); hold on;
                plot(filtered(:,3), filtered(:,6)./filtered(:,3),'b.','MarkerSize',10);
                plot(raw(:,5), raw(:,8)./raw(:,5),'b.','MarkerSize',10);
                ylabel('r (1/s)', 'FontSize', 15);
                xlabel('y (m)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
%                 ylim([-8 2]);
                
                set(gca, 'FontSize', 15); grid on;
                
                plotHandles = opticalExpansionPlot;
                
            else
                assert(length(subplotHandles) == 2, 'Two axes handles are needed to plot the data');
%                 axes(subplotHandles(1));
                plot(filtered(:,3), filtered(:,6),'b.','MarkerSize',10, 'Parent', subplotHandles(1));
                plot(raw(:,5), raw(:,8),'r.','MarkerSize',10, 'Parent', subplotHandles(1));
                legend(subplotHandles(1),{'filtered','raw'},'Location','best');
                                
%                 axes(subplotHandles(2));
                plot(filtered(:,3), filtered(:,6)./filtered(:,3),'b.','MarkerSize',10, 'Parent', subplotHandles(2));
                plot(raw(:,5), raw(:,8)./raw(:,5),'r.','MarkerSize',10, 'Parent', subplotHandles(2));

                plotHandles = [];
            end
        end
        
        
        function plotHandles = plotData(obj)
            if length(obj.rawTrack)~=length(obj.state) % if filtered data doesn't exist, filter it first
                obj.filterRawTracks(20);
            end
            
            trajPlot = figure;
            velPlot = figure;
            accPlot = figure;
            colormap_raw = brewermap(length(obj.rawTrack),'Dark2');
            colormap_filtered = brewermap(length(obj.rawTrack),'Set1');
            t0 = obj.rawTrack(1).rawState(1,3);
            % Plots all raw and filtered track excerpts
            for ct=1:length(obj.rawTrack)
                
                % Plot trajectory
                rawState = obj.rawTrack(ct).rawState;
                filteredState = obj.state(ct).filteredState;
                
                figure(trajPlot);
%                 figHandles(1) = trajPlot;
%                 subplotHandles(1) = subplot(3,1,1); hold on;
                subplot(3,1,1); hold on;
                if ct==1
                    plot(rawState(:,3)-t0,rawState(:,4),'.','MarkerSize',10, 'MarkerFaceColor', colormap_raw(ct,:), 'MarkerEdgeColor', colormap_raw(ct,:), 'DisplayName', 'raw');
                    plot(filteredState(:,1)-t0,filteredState(:,2),'Color', colormap_filtered(ct,:), 'LineWidth', 2, 'DisplayName', 'filtered');
                else
                    plot(rawState(:,3)-t0,rawState(:,4),'.','MarkerSize',10, 'MarkerFaceColor', colormap_raw(ct,:), 'MarkerEdgeColor', colormap_raw(ct,:), 'HandleVisibility', 'off');
                    plot(filteredState(:,1)-t0,filteredState(:,2),'Color', colormap_filtered(ct,:), 'LineWidth', 2, 'HandleVisibility', 'off');
                end
                ylabel('x (m)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                

%             figHandles(end+1) = trajPlot;
%             subplotHandles(end+1) = subplot(3,1,2); hold on;
                subplot(3,1,2); hold on;
                plot(rawState(:,3)-t0,rawState(:,5),'.','MarkerSize',10, 'MarkerFaceColor', colormap_raw(ct,:), 'MarkerEdgeColor', colormap_raw(ct,:));
                plot(filteredState(:,1)-t0,filteredState(:,3),'Color', colormap_filtered(ct,:), 'LineWidth', 2);
                ylabel('y (m)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;

%             figHandles(end+1) = trajPlot;
%             subplotHandles(end+1) = subplot(3,1,3); hold on;
                subplot(3,1,3); hold on;
                plot(rawState(:,3)-t0,rawState(:,6),'.','MarkerSize',10, 'MarkerFaceColor', colormap_raw(ct,:), 'MarkerEdgeColor', colormap_raw(ct,:));
                plot(filteredState(:,1)-t0,filteredState(:,4),'Color', colormap_filtered(ct,:), 'LineWidth', 2);
                ylabel('z (m)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                xlabel('time (s)', 'FontSize', 16);
                
                % Plot velocity
                
                figure(velPlot);
%                 figHandles(end+1) = velPlot;
%                 subplotHandles(end+1) = subplot(3,1,1); hold on;
                subplot(3,1,1); hold on;
                if ct==1
                    plot(rawState(:,3)-t0,rawState(:,7),'.','MarkerSize',10, 'MarkerFaceColor', colormap_raw(ct,:), 'MarkerEdgeColor', colormap_raw(ct,:), 'DisplayName', 'raw');
                    plot(filteredState(:,1)-t0,filteredState(:,5),'Color', colormap_filtered(ct,:), 'LineWidth', 2, 'DisplayName', 'filtered');
                else
                    plot(rawState(:,3)-t0,rawState(:,7),'.','MarkerSize',10, 'MarkerFaceColor', colormap_raw(ct,:), 'MarkerEdgeColor', colormap_raw(ct,:), 'HandleVisibility', 'off');
                    plot(filteredState(:,1)-t0,filteredState(:,5),'Color', colormap_filtered(ct,:), 'LineWidth', 2, 'HandleVisibility', 'off');
                end
                ylabel('V_{gx} (m/s)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
%                 if ct==1
%                     legend({'raw','filtered'});
%                 end
                
%                 figHandles(end+1) = velPlot;
%                 subplotHandles(end+1) = subplot(3,1,2); hold on;
                subplot(3,1,2); hold on;
                plot(rawState(:,3)-t0,rawState(:,8),'.','MarkerSize',10, 'MarkerFaceColor', colormap_raw(ct,:), 'MarkerEdgeColor', colormap_raw(ct,:));
                plot(filteredState(:,1)-t0,filteredState(:,6),'Color', colormap_filtered(ct,:), 'LineWidth', 2);
                ylabel('V_{gy} (m/s)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                
%                 figHandles(end+1) = velPlot;
%                 subplotHandles(end+1) = subplot(3,1,3); hold on;
                subplot(3,1,3); hold on;
                plot(rawState(:,3)-t0,rawState(:,9),'.','MarkerSize',10, 'MarkerFaceColor', colormap_raw(ct,:), 'MarkerEdgeColor', colormap_raw(ct,:));
                plot(filteredState(:,1)-t0,filteredState(:,7),'Color', colormap_filtered(ct,:), 'LineWidth', 2);
                ylabel('V_{gz} (m/s)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                xlabel('time (s)', 'FontSize', 16);
                
                % Plot acceleration
                
                figure(accPlot);
%                 figHandles(end+1) = accPlot;
%                 subplotHandles(end+1) = subplot(3,1,1); hold on;
                subplot(3,1,1); hold on;
                if ct==1
                    plot(filteredState(:,1)-t0,filteredState(:,8),'Color', colormap_filtered(ct,:), 'LineWidth', 2, 'DisplayName', 'filtered');
                else
                    plot(filteredState(:,1)-t0,filteredState(:,8),'Color', colormap_filtered(ct,:), 'LineWidth', 2, 'HandleVisibility', 'off');
                end
%                 plot(filteredState(:,1)-filteredState(1,1),filteredState(:,8),'Color', colormap_filtered(ct,:), 'LineWidth', 2);
                ylabel('a_{gx} (m/s^2)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
%                 if ct==1
%                     legend({'filtered'});
%                 end
                
%                 figHandles(end+1) = accPlot;
%                 subplotHandles(end+1) = subplot(3,1,2); hold on;
                subplot(3,1,2); hold on;
                plot(filteredState(:,1)-t0,filteredState(:,9),'Color', colormap_filtered(ct,:), 'LineWidth', 2);
                ylabel('a_{gy} (m/s^2)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                
%                 figHandles(end+1) = accPlot;
%                 subplotHandles(end+1) = subplot(3,1,3); hold on;
                subplot(3,1,3); hold on;
                plot(filteredState(:,1)-t0,filteredState(:,10),'Color', colormap_filtered(ct,:), 'LineWidth', 2);
                ylabel('a_{gz} (m/s^2)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                xlabel('time (s)', 'FontSize', 16);
            end
            
            figure(trajPlot); subplot(3,1,1);
            legend(gca,'show','Location','best');
            figure(velPlot); subplot(3,1,1);
            legend(gca,'show','Location','best');
            figure(accPlot); subplot(3,1,1);
            legend(gca,'show','Location','best');
            
            plotHandles = [trajPlot; velPlot; accPlot];

%             arrayfun(@(x,y) obj.drawlines(x,y), subplotHandles, figHandles);
            
            obj.displayNext2EachOther([trajPlot; velPlot; accPlot]);
        end
        
        function drawlines(obj, subplotHandle, figHandle)
            figure(figHandle);
            subplot(subplotHandle);
            grid on;
            subplotHandle.GridAlpha = 0.2;
            if ~isempty(obj.state_le)
                vline(obj.state_le(1) - obj.state_le(1),'k');
%                 vline(obj.state_le(1) - obj.state_le(1),'k','LE');
            end
            if ~isempty(obj.state_t1)
                vline(obj.state_t1(1) - obj.state_le(1),'r');
%                 vline(obj.state_t1(1) - obj.state_le(1),'r','T1');
            end
            if ~isempty(obj.state_t2)
                vline(obj.state_t2(1) - obj.state_le(1),'b');
%                 vline(obj.state_t2(1) - obj.state_le(1),'b','T2');
            end
            
        end
        
        function plotHandles = plotData_LDF(obj)
            % Plots trajectory (in landing disc reference frame) with time (le as 0)
            % Also plots Vgy vs y
            % Also plots Vgy/y vs y
            
            if length(obj.rawTrack)~=length(obj.state) % if filtered data doesn't exist, filter it first
                obj.filterRawTracks(20);
            end
            
            if length(obj.state_LDF)~=length(obj.state) % if filtered data doesn't exist, filter it first
                error('First compute the states in Landing Discs Frame...');
            end
            
            trajPlot = figure;
            velPlot = figure;
            accPlot = figure;
            opticalExpansionPlot = figure;
            colormap_filtered = brewermap(length(obj.state_LDF),'Set1');
            t0 = obj.state_LDF(1).filteredState(1,1);
            
            subplotHandles = gobjects(12,1);
            % Plots all track excerpts
            for ct=1:length(obj.state_LDF)
                
                
                filteredState = obj.state_LDF(ct).filteredState;
%                 filteredState(1,1)
                
                % Plot trajectory
                figure(trajPlot);
%                 figHandles(1) = trajPlot;
                subplotHandles(1) = subplot(3,1,1); hold on;
                plot(filteredState(:,1)-t0,filteredState(:,2),'Color', colormap_filtered(ct,:), 'LineWidth', 2);
                ylabel('x (m)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                

%             figHandles(end+1) = trajPlot;
                subplotHandles(2) = subplot(3,1,2); hold on;
%                 subplot(3,1,2); hold on;
                plot(filteredState(:,1)-t0,filteredState(:,3),'Color', colormap_filtered(ct,:), 'LineWidth', 2);
                ylabel('y (m)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;

%             figHandles(end+1) = trajPlot;
                subplotHandles(3) = subplot(3,1,3); hold on;
%                 subplot(3,1,3); hold on;
                plot(filteredState(:,1)-t0,filteredState(:,4),'Color', colormap_filtered(ct,:), 'LineWidth', 2);
                ylabel('z (m)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                xlabel('time (s)', 'FontSize', 16);
                
                linkaxes(subplotHandles(1:3),'x');
                
                % Plot velocity
                
                figure(velPlot);
%                 figHandles(end+1) = velPlot;
                subplotHandles(4) = subplot(3,1,1); hold on;
%                 subplot(3,1,1); hold on;
                if ct==1
                    plot(filteredState(:,1)-t0,filteredState(:,5),'Color', colormap_filtered(ct,:), 'LineWidth', 2, 'DisplayName', 'filtered');
                else
                    plot(filteredState(:,1)-t0,filteredState(:,5),'Color', colormap_filtered(ct,:), 'LineWidth', 2, 'HandleVisibility', 'off');
                end
                ylabel('V_{gx} (m/s)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
%                 if ct==1
%                     legend({'raw','filtered'});
%                 end
                
%                 figHandles(end+1) = velPlot;
                subplotHandles(5) = subplot(3,1,2); hold on;
%                 subplot(3,1,2); hold on;
                plot(filteredState(:,1)-t0,filteredState(:,6),'Color', colormap_filtered(ct,:), 'LineWidth', 2);
                ylabel('V_{gy} (m/s)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                
%                 figHandles(end+1) = velPlot;
                subplotHandles(6) = subplot(3,1,3); hold on;
%                 subplot(3,1,3); hold on;
                plot(filteredState(:,1)-t0,filteredState(:,7),'Color', colormap_filtered(ct,:), 'LineWidth', 2);
                ylabel('V_{gz} (m/s)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                xlabel('time (s)', 'FontSize', 16);
                
                linkaxes(subplotHandles(4:6),'x');
                
                % Plot acceleration
                
                figure(accPlot);
%                 figHandles(end+1) = accPlot;
                subplotHandles(7) = subplot(3,1,1); hold on;
%                 subplot(3,1,1); hold on;
                if ct==1
                    plot(filteredState(:,1)-t0,filteredState(:,8),'Color', colormap_filtered(ct,:), 'LineWidth', 2, 'DisplayName', 'filtered');
                else
                    plot(filteredState(:,1)-t0,filteredState(:,8),'Color', colormap_filtered(ct,:), 'LineWidth', 2, 'HandleVisibility', 'off');
                end
%                 plot(filteredState(:,1)-filteredState(1,1),filteredState(:,8),'Color', colormap_filtered(ct,:), 'LineWidth', 2);
                ylabel('a_{gx} (m/s^2)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
%                 if ct==1
%                     legend({'filtered'});
%                 end
                
%                 figHandles(end+1) = accPlot;
                subplotHandles(8) = subplot(3,1,2); hold on;
%                 subplot(3,1,2); hold on;
                plot(filteredState(:,1)-t0,filteredState(:,9),'Color', colormap_filtered(ct,:), 'LineWidth', 2);
                ylabel('a_{gy} (m/s^2)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                
%                 figHandles(end+1) = accPlot;
                subplotHandles(9) = subplot(3,1,3); hold on;
%                 subplot(3,1,3); hold on;
                plot(filteredState(:,1)-t0,filteredState(:,10),'Color', colormap_filtered(ct,:), 'LineWidth', 2);
                ylabel('a_{gz} (m/s^2)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                xlabel('time (s)', 'FontSize', 16);
                
                linkaxes(subplotHandles(7:9),'x');
                
                
                % Plot optical flow parameters
                % These parameters are plotted only for last 1.5s seconds 
                state_subset = filteredState(filteredState(:,1)-filteredState(end,1)>=-1.5,:);
                
                figure(opticalExpansionPlot);
%                 figHandles(end+1) = opticalExpansionPlot;
                subplotHandles(10) = subplot(3,1,1); hold on;
                plot(state_subset(:,3), state_subset(:,6), '.', ...
                    'MarkerSize', 10, 'MarkerFaceColor', colormap_filtered(ct,:), 'MarkerEdgeColor', colormap_filtered(ct,:));
                ylabel('V_{gy} (m/s)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                
%                 figHandles(end+1) = opticalExpansionPlot;
                subplotHandles(11) = subplot(3,1,2); hold on;
                plot(state_subset(:,3), state_subset(:,6)./state_subset(:,3),'.', ...
                    'MarkerSize', 10, 'MarkerFaceColor', colormap_filtered(ct,:), 'MarkerEdgeColor', colormap_filtered(ct,:));
                ylabel('r (1/s)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                
                subplotHandles(11) = subplot(3,1,2); hold on;
                plot(state_subset(:,3), state_subset(:,6)./state_subset(:,3),'.', ...
                    'MarkerSize', 10, 'MarkerFaceColor', colormap_filtered(ct,:), 'MarkerEdgeColor', colormap_filtered(ct,:));
                ylabel('r (1/s)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                ylim([-10 5]);
                
                subplotHandles(12) = subplot(3,1,3); hold on;
                plot(state_subset(:,3),...
                state_subset(:,9)./state_subset(:,3)-state_subset(:,6).^2./state_subset(:,3).^2, '.', ...
                    'MarkerSize', 10, 'MarkerFaceColor', colormap_filtered(ct,:), 'MarkerEdgeColor', colormap_filtered(ct,:));
                ylabel('r dot (1/s^2)', 'FontSize', 15);
                xlabel('y (m)', 'FontSize', 16);
                set(gca, 'FontSize', 18); grid on;
                linkaxes(subplotHandles(10:12),'x');
                ylim([-100 50]);
                
            end
            
            figure(trajPlot); subplot(3,1,1);
            xlim([0 Inf]);
            figure(velPlot); subplot(3,1,1);
            xlim([0 Inf]);
            figure(accPlot); subplot(3,1,1);
            xlim([0 Inf]);
            
            
            
%             arrayfun(@(x,y,z) obj.drawlines_LDF(x,y,z), subplotHandles, figHandles, xaxis_dim);
            
            obj.displayNext2EachOther([trajPlot; velPlot; accPlot; opticalExpansionPlot]);
            
            plotHandles = [trajPlot; velPlot; accPlot; opticalExpansionPlot];
        end
        
        
        function drawlines_LDF(obj, subplotHandle, figHandle, xaxis_dim)
            
            figure(figHandle);
            subplot(subplotHandle);
            grid on;
%             subplotHandle.GridAlpha = 0.2;
            if ~isempty(obj.state_le_LDF) && strcmpi(xaxis_dim, 'time')
                vline(obj.state_le_LDF(1) - obj.state_le_LDF(1),'k');
%                 vline(obj.state_le_LDF(1) - obj.state_le_LDF(1),'k','LE');
            elseif ~isempty(obj.state_le_LDF) && strcmpi(xaxis_dim, 'distance')
                vline(obj.state_le_LDF(3),'k');
            end
            if ~isempty(obj.state_t1) && strcmpi(xaxis_dim, 'time')
                vline(obj.state_t1_LDF(1) - obj.state_le_LDF(1),'r');
%                 vline(obj.state_t1_LDF(1) - obj.state_le_LDF(1),'r','T1');
            elseif ~isempty(obj.state_t1_LDF) && strcmpi(xaxis_dim, 'distance')
                vline(obj.state_t1_LDF(3),'r');
            end
            if ~isempty(obj.state_t2) && strcmpi(xaxis_dim, 'time')
                vline(obj.state_t2_LDF(1) - obj.state_le_LDF(1),'b');
%                 vline(obj.state_t2_LDF(1) - obj.state_le_LDF(1),'b','T2');
            elseif ~isempty(obj.state_t2_LDF) && strcmpi(xaxis_dim, 'distance')
                vline(obj.state_t2_LDF(3),'b');
            end
            
        end
        
        function displayNext2EachOther(obj, plotHandles)
            pos_screen = get(0,'ScreenSize');
            for ct=1:length(plotHandles)
                if ct==1
                    plotHandles(ct).Position = [1 pos_screen(4)/length(plotHandles) plotHandles(ct).Position(3:4) ];
                else
                    plotHandles(ct).Position = [(ct-1)*pos_screen(3)/length(plotHandles) pos_screen(4)/length(plotHandles) plotHandles(ct).Position(3:4) ];
                end

            end
        end
        
        function drawVvsy(obj, axesHandle)
            % Draws V_{gy} vs y plot on plotHandle
            % axesHandle - handle of axes from GUI figure
            % Used for plotting in GUI for track extraction for sysID
            state_subset = obj.state_LDF(obj.time>=-1.5,:);
            plot(axesHandle, state_subset(:,3),state_subset(:,6),'m.','MarkerSize',10);
            axesHandle.XLim = [0 0.2];
        end
        function drawVbyyvsy(obj, axesHandle)
            % Draws V_{gy}/y vs y plot on plotHandle
            % axesHandle - handle of axes from GUI figure
            % Used for plotting in GUI for track extraction for sysID
            state_subset = obj.state_LDF(obj.time>=-1.5,:);
            plot(axesHandle, state_subset(:,3),state_subset(:,6)./state_subset(:,3),'m.','MarkerSize',10);
            axesHandle.XLim = [0 0.2];
        end
        
        function drawOpticalExpansion(obj, plotHandle)
            state_subset = obj.state_LDF(obj.time>=-1.5,:);
            figure(plotHandle);
            p1 = subplot(2,1,1); hold on;
            plot(state_subset(:,3),state_subset(:,6),'m.','MarkerSize',10);
            ylabel('V_{gy} (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;
            
            p2 = subplot(2,1,2); hold on;
            plot(state_subset(:,3),state_subset(:,6)./state_subset(:,3),'m.','MarkerSize',10);
            ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
            xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;
            
            linkaxes([p1, p2],'x');
        end
        
        function fftHandle = plotFFT(obj, time_duration)
            % time_duration - # of seconds before T1 for which FFT is computed
            
            time_duration = 2.00188159942627; % corresponds to 350 samples
            
            fs = 1/obj.dt;
%             L = round(time_duration/obj.dt);
            
            response = obj.state((obj.state(:,1)-obj.state_t1(1)>=-(time_duration+1e-5)) & ...
                          (obj.state(:,1)-obj.state_t1(1)<=1e-5), [3 6 9]);
            L = size(response, 1);
            
            if rem(L,2)~=0
                response(1,:) = [];
                L = L-1;
            end
            
            
            response_fft = fft(response);
            P2 = abs(response_fft/L);
            P1 = P2(1:L/2+1,:);
            P1(2:end-1,:) = 2*P1(2:end-1,:);
            
            f = fs*(0:(L/2))/L;
            fftHandle = figure;
            subplot(3,1,1);
            plot(f,P1(:,1));
            xlabel('f (Hz)', 'FontSize', 16);
            ylabel('|P_y(f)|', 'FontSize', 16);
            
            subplot(3,1,2);
            plot(f,P1(:,2));
            xlabel('f (Hz)', 'FontSize', 16);
            ylabel('|P_{Vy}(f)|', 'FontSize', 16);
            
            subplot(3,1,3);
            plot(f,P1(:,1));
            xlabel('f (Hz)', 'FontSize', 16);
            ylabel('|P_{ay}(f)|', 'FontSize', 16);
            
        end
        
        function plotHandles = plotData_Distance(obj)
            % Plot optical flow parameters
            % These parameters are plotted only for 1.5s seconds before LE
            
            str = '#375E98';
            color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

            state_subset = obj.state_LDF(obj.time>=-1.5,:);
            
            opticalExpansionPlot = figure;
            figure(opticalExpansionPlot);
            figHandles(1) = opticalExpansionPlot;
            subplotHandles(1) = subplot(3,1,1); hold on;
            xaxis_dim{1} = 'distance';
%             plot(obj.state_LDF(:,3),obj.state_LDF(:,6),'m.','MarkerSize',10);
            plot(state_subset(:,3),state_subset(:,6),'b.','MarkerSize',10);
            ylabel('V_{gy} (m/s)', 'FontSize', 15);
%             ylabel('Velocity towards the platform (m/s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            
            
            figHandles(end+1) = opticalExpansionPlot;
            subplotHandles(end+1) = subplot(3,1,2); hold on;
            xaxis_dim{end+1} = 'distance';
%             plot(obj.state_LDF(:,3),obj.state_LDF(:,6)./obj.state_LDF(:,3),'m.','MarkerSize',10);
            plot(state_subset(:,3),state_subset(:,6)./state_subset(:,3),'b.','MarkerSize',10);
            ylabel('r (1/s)', 'FontSize', 15);
%             ylabel('Optical rate of expansion (1/s)', 'FontSize', 15);
%             xlabel('Distance from the platform (m)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            ylim([-10 1]);
            
            M = movvar(state_subset(:,6)./state_subset(:,3), [174 1]);
            figHandles(end+1) = opticalExpansionPlot;
            subplotHandles(end+1) = subplot(3,1,3); hold on;
            xaxis_dim{end+1} = 'distance';
%             plot(obj.state_LDF(:,3),obj.state_LDF(:,6)./obj.state_LDF(:,3),'m.','MarkerSize',10);
            plot(state_subset(:,3),M,'b.','MarkerSize',10);
            ylabel('Var. of r', 'FontSize', 15);
            xlabel('Distance from the platform (m)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            ylim([0 4]);
            
            linkaxes(subplotHandles(1:3),'x');
            xlim([0 0.2]);
            
            arrayfun(@(x,y,z) obj.drawlines_Poster(x,y,z), subplotHandles, figHandles, xaxis_dim);
                        
            plotHandles = [opticalExpansionPlot];
        end
        
        
        
        function drawlines_Poster(obj, subplotHandle, figHandle, xaxis_dim)
            
            figure(figHandle);
            subplot(subplotHandle);
%             grid on;
%             subplotHandle.GridAlpha = 0.2;
%             if ~isempty(obj.state_le_LDF) && strcmpi(xaxis_dim, 'time')
%                 vline(obj.state_le_LDF(1) - obj.state_le_LDF(1),'k');
% %                 vline(obj.state_le_LDF(1) - obj.state_le_LDF(1),'k','LE');
%             elseif ~isempty(obj.state_le_LDF) && strcmpi(xaxis_dim, 'distance')
%                 vline(obj.state_le_LDF(3),'k');
%             end
            if ~isempty(obj.state_t1) && strcmpi(xaxis_dim, 'time')
                vline(obj.state_t1_LDF(1) - obj.state_le_LDF(1),'r');
%                 vline(obj.state_t1_LDF(1) - obj.state_le_LDF(1),'r','T1');
            elseif ~isempty(obj.state_t1_LDF) && strcmpi(xaxis_dim, 'distance')
                vline(obj.state_t1_LDF(3),'r');
            end
%             if ~isempty(obj.state_t2) && strcmpi(xaxis_dim, 'time')
%                 vline(obj.state_t2_LDF(1) - obj.state_le_LDF(1),'b');
% %                 vline(obj.state_t2_LDF(1) - obj.state_le_LDF(1),'b','T2');
%             elseif ~isempty(obj.state_t2_LDF) && strcmpi(xaxis_dim, 'distance')
%                 vline(obj.state_t2_LDF(3),'b');
%             end
            
        end
        
        function filterTrack_old(obj)
%             % To filter the raw tracks obtained from FLYDRA
%             for ct_id=1:obj.nObjects
%                 obj_id = obj.obj_ids(ct_id);
%                 
%                 track = sortrows(obj.rawTrack(obj.rawTrack(:,1)==obj_id,1:9), [3 2]);
%                 
%                 % Averaging multiple values at the same timestamp
%                 % (delta_t<1e-5)
%                 diff_t = [1/175; diff(track(:,3))];
%                 unique_t_indices = [find(abs(diff_t)>1e-5); size(track,1)+1];
%                 unique_track = zeros(numel(unique_t_indices)-1,9);
%                 for i=1:length(unique_t_indices)-1
%                     unique_track(i,:) = mean(track(unique_t_indices(i):unique_t_indices(i+1)-1,:),1);
%                 end
%                 
%                 obj.state = zeros(length(unique_track),10);
%                 obj.state(:,1:7) = unique_track(:,3:9);
% %                 obj.state(:,1:4) = unique_track(:,3:6);
% %                 obj.state(:,5:7) = diffxy(obj.state(:,1), obj.state(:,2:4));
%                 obj.state(:,8:10) = diffxy(obj.state(:,1), obj.state(:,5:7));
%                 % Filtering and computing derivatives based on the filtered
%                 % signal
%                 % Signal is divided into subparts if more than 4
%                 % consecutive values are missing. Filtering is applied to
%                 % each of those parts separately
%                 obj.state(:,5:7) = diffxy(obj.state(:,1), obj.state(:,2:4));
%                 obj.state(:,8:10) = diffxy(obj.state(:,1), obj.state(:,5:7));
%                 
%                 diff_t = diff(unique_track(:,3));
%                 parts = [1; find(abs(diff_t)>4*1/175); size(unique_track,1)];
%                 
%                 for i=1:length(parts)-1
%                     time = unique_track(parts(i):parts(i+1),3);
%                     mean_dt = mean(diff(time));
%                     
%                     
% %                     mean_dt = mean(diff_t(parts(i):))
%                 end
%                 
%                 
%                 [track(:,3) diff_t abs(diff_t)>1e-5]
%                 
%                 
%                 
%                 
%                 [num,den]=butter(5, fc/(fs/2),'low'); % 5-th order Butterworth filter %XXXX CHECK ORDER OF THE FILTER
%             end
%             
%             
        end
        
    end
    
    methods(Static)
        function state_LDF = convert_to_landing_disc_reference_frame1(state, origin, landing_side)
            % state - N X 10 (or N X 7 in case of rawData) matrix similar to Landingtrack.state
            % origin - 1 X 3 vector containing origin of the landing
            % disc reference frame
            % landing_side - 'Hive' or 'Feeder'
            
            % returns state (N X 10) in landing disc reference frame (state_LDF)
                
            if isempty(state)
                state_LDF = state;
                
            elseif strcmpi(landing_side,'Hive')
                % This reference frame is left handed with x in the
                % downstream direction, y pointing towards the feeder and z
                % pointing up. The origin is at the center of the landing
                % disc
                
                state_LDF = state;
                state_LDF(:,2:4) = state(:,2:4) - origin;
                state_LDF(:, [3 6 9]) = -state_LDF(:, [3 6 9]);

            elseif strcmpi(landing_side,'Feeder')
                % This reference frame is right handed with x in the
                % downstream direction, y pointing towards the hive and z
                % pointing up. The origin is at the center of the landing
                % disc
                
                state_LDF = state;
                state_LDF(:,2:4) = state(:,2:4) - origin;
                
            else
                error('Unknown landing disc found :/. How is that possible? MAGIC??');
            end
            
        end
        
        function state_LDF = convert_to_landing_disc_reference_frame2(state, origin, landing_side)
            % state - N X 10 matrix similar to Landingtrack.state
            % origin - 1 X 3 vector containing origin of the landing
            % disc reference frame
            % landing_side - 'Hive' or 'Feeder'
            
            % returns state (N X 10) in landing disc reference frame (state_LDF)
                
            if isempty(state)
                state_LDF = state;
                
            elseif strcmpi(landing_side,'Hive')
                % This reference frame is right handed with x in the
                % downstream direction, y pointing in the direction opposite to the feeder and z
                % pointing up. The origin is at the center of the landing
                % disc
                
                state_LDF = state;
                state_LDF(:,2:4) = state(:,2:4) - origin;

            elseif strcmpi(landing_side,'Feeder')
                % This reference frame is left handed with x in the
                % downstream direction, y pointing in the direction opposite to the hive and z
                % pointing up. The origin is at the center of the landing
                % disc
                
                state_LDF = state;
                state_LDF(:,2:4) = state(:,2:4) - origin;
                if size(state,2) <= 7
                    state_LDF(:, [3 6]) = -state_LDF(:, [3 6]);
                else
                    state_LDF(:, [3 6 9]) = -state_LDF(:, [3 6 9]);
                end
                
            else
                error('Unknown landing disc found :/. How is that possible? MAGIC??');
            end
            
        end
        
        function plotHandles = plotDataLDF_Time_rdot(state_LDF, subplotHandles)
            % Plot optical flow parameters
            % These parameters are plotted with time
            
            % state_LDF - instance of filteredState_BlindLandingTrack
            % subplotHandles - array of axes handles to plot data in
           

            filteredState = state_LDF.filteredState;
            state_subset = filteredState; %filteredState(filteredState(:,1)-filteredState(end,1)>=-1.3,:);
            t0 = filteredState(end,1);
            
            if nargin == 1
                opticalExpansionPlot = figure;
                figure(opticalExpansionPlot);
                subplotHandles(1) = subplot(5,1,1); hold on;
                plot(state_subset(:,1)-t0, state_subset(:,6),'b.','MarkerSize',10);
                ylabel('V_{gy} (m/s)', 'FontSize', 15);
    %             ylabel('Velocity towards the platform (m/s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;

                subplotHandles(end+1) = subplot(5,1,2); hold on;
                plot(state_subset(:,1)-t0, state_subset(:,3),'b.','MarkerSize',10);
                ylabel('y (m)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
    %             ylim([-10 1]);

                subplotHandles(end+1) = subplot(5,1,3); hold on;
                plot(state_subset(:,1)-t0, state_subset(:,6)./state_subset(:,3),'b.','MarkerSize',10);
                ylabel('r (1/s)', 'FontSize', 15);
    %             ylabel('Optical rate of expansion (1/s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                ylim([-10 2]);

                subplotHandles(end+1) = subplot(5,1,4); hold on;
                rdot = diffxy(state_subset(:,1)-state_subset(1,1), state_subset(:,6)./state_subset(:,3));
                plot(state_subset(:,1)-t0,rdot,'b.','MarkerSize',10); hold on;
                plot(state_subset(:,1)-t0,...
                    state_subset(:,9)./state_subset(:,3)-state_subset(:,6).^2./state_subset(:,3).^2,'r.','MarkerSize',10);
                ylabel('r dot (1/s^2)', 'FontSize', 15);
    %             xlabel('Time (s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                ylim([-50 30]);

                subplotHandles(end+1) = subplot(5,1,5); hold on;
                plot(state_subset(:,1)-t0,...
                    state_subset(:,9),'b.','MarkerSize',10);
                ylabel('a_y (m/s^2)', 'FontSize', 15);
                xlabel('Time (s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                
                linkaxes(subplotHandles, 'x');
                plotHandles = [opticalExpansionPlot];
            elseif nargin == 2
                assert(length(subplotHandles) == 5, 'Five axes handles are needed to plot the data');
                axes(subplotHandles(1));
                plot(state_subset(:,1)-t0, state_subset(:,6),'b.','MarkerSize',10);
                
                axes(subplotHandles(2));
                plot(state_subset(:,1)-t0, state_subset(:,3),'b.','MarkerSize',10);
                
                axes(subplotHandles(3));
                plot(state_subset(:,1)-t0, state_subset(:,6)./state_subset(:,3),'b.','MarkerSize',10);
                
                axes(subplotHandles(4));
                plot(state_subset(:,1)-t0,...
                    state_subset(:,9)./state_subset(:,3)-state_subset(:,6).^2./state_subset(:,3).^2,'b.','MarkerSize',10);
                
                axes(subplotHandles(5));
                plot(state_subset(:,1)-t0,...
                    state_subset(:,9),'b.','MarkerSize',10);
                
                plotHandles = [];
            end
        end
        
        function plotHandles = plotDataLDF_Distance_rdot(state_LDF, subplotHandles)
            % Plot optical flow parameters
            % These parameters are plotted with distance from the platform
            
            % state_LDF - state in landing disc frame (as an instance of
            % filteredState_BlindLandingTrack)
            % subplotHandles - array of axes handles to plot data in
            
%             str = '#375E98';
%             color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

            filteredState = state_LDF.filteredState;
            state_subset = filteredState; %filteredState(filteredState(:,1)-filteredState(end,1)>=-1.3,:);
            
            if nargin == 1
                opticalExpansionPlot = figure;
                figure(opticalExpansionPlot);
                subplotHandles(1) = subplot(4,1,1); hold on;
                plot(state_subset(:,3), state_subset(:,6),'b.','MarkerSize',10);
                ylabel('V_{gy} (m/s)', 'FontSize', 15);
    %             ylabel('Velocity towards the platform (m/s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;

                subplotHandles(end+1) = subplot(4,1,2); hold on;
                plot(state_subset(:,3), state_subset(:,6)./state_subset(:,3),'b.','MarkerSize',10);
                ylabel('r (1/s)', 'FontSize', 15);
    %             ylabel('Optical rate of expansion (1/s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                ylim([-8 2]);

                subplotHandles(end+1) = subplot(4,1,3); hold on;
                plot(state_subset(:,3),...
                    state_subset(:,9)./state_subset(:,3)-state_subset(:,6).^2./state_subset(:,3).^2,'b.','MarkerSize',10);
                ylabel('r dot (1/s^2)', 'FontSize', 15);
    %             xlabel('Time (s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                ylim([-50 30]);

                subplotHandles(end+1) = subplot(4,1,4); hold on;
                plot(state_subset(:,3),...
                    state_subset(:,9),'b.','MarkerSize',10);
                ylabel('a_y (m/s^2)', 'FontSize', 15);
                xlabel('Distance from the platform (m)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                
                plotHandles = opticalExpansionPlot;
                
            elseif nargin == 2
                assert(length(subplotHandles) == 4, 'Four axes handles are needed to plot the data');
                axes(subplotHandles(1));
                plot(state_subset(:,3), state_subset(:,6),'b.','MarkerSize',10);
                                
                axes(subplotHandles(2));
                plot(state_subset(:,3), state_subset(:,6)./state_subset(:,3),'b.','MarkerSize',10);
                
                axes(subplotHandles(3));
                plot(state_subset(:,3),...
                    state_subset(:,9)./state_subset(:,3)-state_subset(:,6).^2./state_subset(:,3).^2,'b.','MarkerSize',10);
                
                axes(subplotHandles(4));
                plot(state_subset(:,3),...
                    state_subset(:,9),'b.','MarkerSize',10);
                
                plotHandles = [];
            end
            
            
            
            
%             M = movvar(state_subset(:,6)./state_subset(:,3), [9 0]);
%             figHandles(end+1) = opticalExpansionPlot;
%             subplotHandles(end+1) = subplot(4,1,4); hold on;
%             xaxis_dim{end+1} = 'time';
%             plot(state_subset(:,1)-obj.state_le_LDF(1,1),M,'b.','MarkerSize',10);
%             ylabel('Var. of r', 'FontSize', 15);
%             xlabel('Time to leg extension (s)', 'FontSize', 15);
%             set(gca, 'FontSize', 15); grid on;
%             ylim([0 4]);
            
            linkaxes(subplotHandles(1:end),'x');
%             xlim([-1.5 0]);
                        
            
        end
        
        function plotHandles = plotDataLDF_Time_forArticle(state_LDF, t_start)
            % Plot optical flow parameters
            % These parameters are plotted with time
            
            % state_LDF - instance of filteredState_BlindLandingTrack
            % t_start - 
           

            filteredState = state_LDF.filteredState;
            state_subset = filteredState; % filteredState(filteredState(:,1)-filteredState(end,1)>=-1.3,:);
            t0 = filteredState(find(filteredState(:,1)>t_start,1),1);
            tend = filteredState(end,1);
            
            if nargin == 2
                opticalExpansionPlot = figure;
                figure(opticalExpansionPlot);
                
                subplotHandles(1) = subplot(3,1,1); hold on;
%                 plot(state_subset(:,1)-t0, -1*state_subset(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                plot(state_subset(:,1)-tend, -1*state_subset(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                ylabel('y (m)', 'FontSize', 16);
                set(gca, 'FontSize', 15); grid on;
                
                subplotHandles(2) = subplot(3,1,2); hold on;
%                 plot(state_subset(:,1)-t0, state_subset(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                plot(state_subset(:,1)-tend, state_subset(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                ylabel('V (m/s)', 'FontSize', 16);
    %             ylabel('Optical rate of expansion (1/s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                
                subplotHandles(3) = subplot(3,1,3); hold on;
%                 plot(state_subset(:,1)-t0, -1*state_subset(:,6)./state_subset(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                plot(state_subset(:,1)-tend, -1*state_subset(:,6)./state_subset(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                ylabel('r (1/s)', 'FontSize', 16);
                xlabel('time (s)', 'FontSize', 16);
                set(gca, 'FontSize', 15); grid on;

                plotHandles = [opticalExpansionPlot];
                linkaxes(subplotHandles(1:3),'x');
            end
        end
        
        function plotHandles = plotDataLDF_Time2_forArticle(state_LDF, t_start)
            % Plot optical flow parameters
            % These parameters are plotted with time
            
            % state_LDF - instance of filteredState_BlindLandingTrack
            % t_start - 
           

            filteredState = state_LDF.filteredState;
            state_subset = filteredState; % filteredState(filteredState(:,1)-filteredState(end,1)>=-1.3,:);
            t0 = filteredState(find(filteredState(:,1)>t_start,1),1);
            tend = filteredState(end,1);
            
            if nargin == 2
                opticalExpansionPlot = figure;
                figure(opticalExpansionPlot);
                
                subplotHandles(1) = subplot(3,1,1); hold on;
%                 plot(state_subset(:,1)-t0, -1*state_subset(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                plot(state_subset(:,1)-tend, -1*state_subset(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                ylabel('y (m)', 'FontSize', 16);
                set(gca, 'FontSize', 15); grid on;
                
                subplotHandles(2) = subplot(3,1,2); hold on;
%                 plot(state_subset(:,1)-t0, state_subset(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                plot(state_subset(:,1)-tend, state_subset(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                ylabel('V (m/s)', 'FontSize', 16);
    %             ylabel('Optical rate of expansion (1/s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                
                subplotHandles(3) = subplot(3,1,3); hold on;
%                 plot(state_subset(:,1)-t0, -1*state_subset(:,6)./state_subset(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                plot(state_subset(:,1)-tend, -1*state_subset(:,3)./state_subset(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                ylabel('\tau (s)', 'FontSize', 16);
                xlabel('time (s)', 'FontSize', 16);
                set(gca, 'FontSize', 15); grid on;

                plotHandles = [opticalExpansionPlot];
                linkaxes(subplotHandles(1:3),'x');
            end
        end
        
        function plotHandles = plotDataLDF_Time3_forArticle(state_LDF, t_start)
            % Plot optical flow parameters
            % These parameters are plotted with time
            
            % state_LDF - instance of filteredState_BlindLandingTrack
            % t_start - 
           

            filteredState = state_LDF.filteredState;
            state_subset = filteredState; % filteredState(filteredState(:,1)-filteredState(end,1)>=-1.3,:);
%             t0 = filteredState(find(filteredState(:,1)>t_start,1),1);
            state_subset = filteredState(filteredState(:,1)>=t_start,:);
            
            tend = filteredState(end,1);
            
            if nargin == 2
                opticalExpansionPlot = figure;
                figure(opticalExpansionPlot);
                
                subplotHandles(1) = subplot(4,1,1); hold on; grid off;
%                 plot(state_subset(:,1)-t0, -1*state_subset(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                plot(state_subset(:,1)-tend, -1*state_subset(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                ylabel('y (m)', 'FontSize', 16);
                set(gca, 'FontSize', 15); grid on;
                
                subplotHandles(2) = subplot(4,1,2); hold on; grid off;
%                 plot(state_subset(:,1)-t0, state_subset(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                plot(state_subset(:,1)-tend, state_subset(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                ylabel('V (ms-1)', 'FontSize', 16);
    %             ylabel('Optical rate of expansion (1/s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                
                subplotHandles(3) = subplot(4,1,3); hold on; grid off;
                plot(state_subset(:,1)-tend, state_subset(:,9),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                ylabel('a (ms-2)', 'FontSize', 16);
%                 xlabel('time (s)', 'FontSize', 16);
                set(gca, 'FontSize', 15); grid on;
                yline(0,'--');
                
                subplotHandles(4) = subplot(4,1,4); hold on; grid off;
%                 plot(state_subset(:,1)-t0, -1*state_subset(:,6)./state_subset(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                plot(state_subset(:,1)-tend, -1*state_subset(:,6)./state_subset(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
                ylabel('r (s-1)', 'FontSize', 16);
                xlabel('time (s)', 'FontSize', 16);
                set(gca, 'FontSize', 15); grid on;

                
                
                
                plotHandles = [opticalExpansionPlot];
                linkaxes(subplotHandles(1:4),'x');
                xlim([t_start-tend 0]);

            end
        end
        
        function plotHandles = plotDataLDF_Time(state_LDF, subplotHandles)
            % Plot optical flow parameters
            % These parameters are plotted with time
            
            % state_LDF - instance of filteredState_BlindLandingTrack
            % subplotHandles - array of axes handles to plot data in
           

            filteredState = state_LDF.filteredState;
            state_subset = filteredState; %filteredState(filteredState(:,1)-filteredState(end,1)>=-1.3,:);
            t0 = filteredState(end,1);
            
            if nargin == 1
                opticalExpansionPlot = figure;
                figure(opticalExpansionPlot);
                
                subplotHandles(1) = subplot(3,1,1); hold on;
                plot(state_subset(:,1)-t0, state_subset(:,6)./state_subset(:,3),'b.','MarkerSize',10);
                ylabel('V_{gy} / y (1/s)', 'FontSize', 15);
    %             ylabel('Optical rate of expansion (1/s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                ylim([-8 2]);

                subplotHandles(2) = subplot(3,1,2); hold on;
                plot(state_subset(:,1)-t0,...
                    state_subset(:,9)./state_subset(:,6),'b.','MarkerSize',10);
                ylabel('a_y / V_{gy} (1/s)', 'FontSize', 15);
    %             xlabel('Time (s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                ylim([-60 60]);

                subplotHandles(3) = subplot(3,1,3); hold on;
                plot(state_subset(:,1)-t0,...
                    state_subset(:,9),'b.','MarkerSize',10);
                ylabel('a_y (m/s^2)', 'FontSize', 15);
                xlabel('Time (s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                
                plotHandles = [opticalExpansionPlot];
                
            elseif nargin == 2
                
                assert(length(subplotHandles) == 3, 'Three axes handles are needed to plot the data');
                
%                 set(0,'CurrentFigure',subplotHandles(1).Parent);
%                 axes(subplotHandles(1));
                plot(state_subset(:,1)-t0, state_subset(:,6)./state_subset(:,3),'b.','MarkerSize',10, 'Parent', subplotHandles(1));
%                 ylim([-8 2]);
                
%                 axes(subplotHandles(2));
                plot(state_subset(:,1)-t0, state_subset(:,9)./state_subset(:,6),'b.','MarkerSize',10, 'Parent', subplotHandles(2));
%                 ylim([-60 60]);
                
%                 axes(subplotHandles(3));
                plot(state_subset(:,1)-t0, state_subset(:,9),'b.','MarkerSize',10, 'Parent', subplotHandles(3));
                
                plotHandles = [];
            end
        end
        
        
        
        function plotHandles = plotDataLDF_abyV(state_LDF, figureHandles)
            % Plot optical flow parameters
            % a vs V
            
            % state_LDF - instance of filteredState_BlindLandingTrack
            % subplotHandles - array of axes handles to plot data in
           

            filteredState = state_LDF.filteredState;
            state_subset = filteredState; %filteredState(filteredState(:,1)-filteredState(end,1)>=-1.3,:);
            
            if nargin == 1
                opticalExpansionPlot = figure;
                figure(opticalExpansionPlot);

                plot(state_subset(:,6), state_subset(:,9), '-b.', 'MarkerSize', 10);
                xlabel('V_{gy} (m/s)', 'FontSize', 15);
                ylabel('a_y (1/s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
%                 ylim([-8 2]);
                
                plotHandles = [opticalExpansionPlot];
                
            elseif nargin == 2
                assert(length(figureHandles) == 1, 'one figure handle is needed to plot the data');
%                 figure(figureHandles(1));
                plot(state_subset(:,6), state_subset(:,9), '-b.', 'MarkerSize', 10, 'Parent', figureHandles(1).Children(1));
                ylim([-30 30]);
                
                plotHandles = [];
            end
        end
        
        function plotHandles = plotDataLDF_Distance(state_LDF, subplotHandles)
            % Plot optical flow parameters
            % These parameters are plotted with distance from the platform
            
            % state_LDF - state in landing disc frame (as an instance of
            % filteredState_BlindLandingTrack)
            % subplotHandles - array of axes handles to plot data in
            
%             str = '#375E98';
%             color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

            filteredState = state_LDF.filteredState;
            state_subset = filteredState; %filteredState(filteredState(:,1)-filteredState(end,1)>=-1.3,:);
            
            if nargin == 1 % only state_LDF is given
                opticalExpansionPlot = figure;
                figure(opticalExpansionPlot);
                subplot(4,1,1); hold on;
                plot(state_subset(:,3), state_subset(:,6),'b.','MarkerSize',10);
                ylabel('V_{gy} (m/s)', 'FontSize', 15);
    %             ylabel('Velocity towards the platform (m/s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;

                subplot(4,1,2); hold on;
                plot(state_subset(:,3), state_subset(:,6)./state_subset(:,3),'b.','MarkerSize',10);
                ylabel('V_{gy} / y (1/s)', 'FontSize', 15);
    %             ylabel('Optical rate of expansion (1/s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                ylim([-8 2]);

                subplot(4,1,3); hold on;
                plot(state_subset(:,3),...
                     state_subset(:,9)./state_subset(:,6),'b.','MarkerSize',10);
                ylabel('a / V_{gy} (1/s)', 'FontSize', 15);
    %             xlabel('Time (s)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                ylim([-50 50]);

                subplot(4,1,4); hold on;
                plot(state_subset(:,3),...
                    state_subset(:,9),'b.','MarkerSize',10);
                ylabel('a_y (m/s^2)', 'FontSize', 15);
                xlabel('y (m)', 'FontSize', 15);
%                 xlabel('Distance from the platform (m)', 'FontSize', 15);
                set(gca, 'FontSize', 15); grid on;
                
                plotHandles = opticalExpansionPlot;
                
            elseif nargin == 2
                assert(length(subplotHandles) == 4, 'Four axes handles are needed to plot the data');
%                 axes(subplotHandles(1));
                plot(state_subset(:,3), state_subset(:,6),'b.','MarkerSize',10, 'Parent', subplotHandles(1));
                                
%                 axes(subplotHandles(2));
                plot(state_subset(:,3), state_subset(:,6)./state_subset(:,3),'b.','MarkerSize',10, 'Parent', subplotHandles(2));
                
%                 axes(subplotHandles(3));
                plot(state_subset(:,3),...
                     state_subset(:,9)./state_subset(:,6),'b.','MarkerSize',10, 'Parent', subplotHandles(3));
                
%                 axes(subplotHandles(4));
                plot(state_subset(:,3),...
                     state_subset(:,9),'b.','MarkerSize',10, 'Parent', subplotHandles(4));
                
                plotHandles = [];
            end
        end
        
        function [datenum, pattern, patternnum, beeID, setID, flightID] = extractInfo(basedir, filename)
            % Function to extract info from the filename in basedir for
            % honeybee data. This code is written based on the readme files
            % in honeybee data folders.
            if strcmpi(basedir, 'Big and Small spiral')
                day_of_exp = datetime(filename(1:6), 'InputFormat', 'ddMMyy');
                datenum = str2num(datestr(day_of_exp, 'yyyymmdd'));
                
                if strcmpi(filename(12:13),'bs')
                    pattern = 'Big static spiral';
%                     patternnum = 1;
                    patternnum = 8; % same as static 4-arm spiral
                elseif strcmpi(filename(12:13),'ss')
                    pattern = 'Small static spiral';
                    patternnum = 2;
                else
                    warning('Unknown Pattern encountered');
                    keyboard;
                end
                
                setID = str2num(filename(15));
                
                flightID = str2num(filename(17:18));
                
                beeID = filename(20:21);
                
            elseif strcmpi(basedir, 'Check') || strcmpi(basedir, 'Grey pattern') || strcmpi(basedir, 'Sector pattern') || strcmpi(basedir, 'Concentric circles')
                
                if strcmpi(filename(1:2),'ch')
                    pattern = 'Checkerboard';
                    patternnum = 3;
                    datenum = 99999999;
                elseif strcmpi(filename(1:2),'gs')
                    pattern = 'Grey pattern';
                    patternnum = 4;
                    datenum = 99999999;
                elseif strcmpi(filename(1:2),'ss')
                    pattern = 'Spoke pattern';
                    patternnum = 11;
                    datenum = 99999999;
                elseif strcmpi(filename(1:2),'cs')
                    pattern = 'Concentric circles';
                    patternnum = 18;
                    datenum = 99999999;
                else
                    warning('Unknown Pattern encountered');
                    keyboard;
                end
                
                split_filename = split(filename, ' ');
                
                setID = str2num(split_filename{2});
                
                flightID = str2num(split_filename{3});
                
                beeID = split_filename{4}(1:2);
                
            elseif strcmpi(basedir, 'Horizontal stripes') || strcmpi(basedir, 'Random Julesz pattern') ...
                    || strcmpi(basedir, 'Vertical stripes')
                day_of_exp = datetime(filename(1:6), 'InputFormat', 'ddMMyy');
                datenum = str2num(datestr(day_of_exp, 'yyyymmdd'));
                
                if strcmpi(basedir, 'Horizontal stripes')
                    pattern = 'Horizontal stripes';
                    patternnum = 5;
                elseif strcmpi(basedir, 'Random Julesz pattern')
                    pattern = 'Random Julesz';
                    patternnum = 6;
                elseif strcmpi(basedir, 'Vertical stripes')
                    pattern = 'Vertical stripes';
                    patternnum = 10;
                else
                    warning('Unknown Pattern encountered');
                    keyboard;
                end
                    
                split_filename = split(filename, ' ');
                
                setID = str2num(split_filename{2});
                
                flightID = str2num(split_filename{3});
                
                beeID = split_filename{4}(1:2);
                
            elseif strcmpi(basedir, 'Variation in number of spiral arms')
                datenum = 99999999;
                
                if strcmpi(filename(1:2), 'S3')
                    pattern = 'Static 3-arm spiral';
                    patternnum = 7;
                elseif strcmpi(filename(1:2), 'S4')
                    pattern = 'Static 4-arm spiral';
                    patternnum = 8;
                elseif strcmpi(filename(1:2), 'S6')
                    pattern = 'Static 6-arm spiral';
                    patternnum = 9;
                elseif strcmpi(filename(1:2), 'E3')
                    pattern = 'Expanding 3-arm spiral';
                    patternnum = 12;
                elseif strcmpi(filename(1:2), 'E4')
                    pattern = 'Expanding 4-arm spiral';
                    patternnum = 13;
                elseif strcmpi(filename(1:2), 'E6')
                    pattern = 'Expanding 6-arm spiral';
                    patternnum = 14;
                elseif strcmpi(filename(1:2), 'C3')
                    pattern = 'Contracting 3-arm spiral';
                    patternnum = 15;
                elseif strcmpi(filename(1:2), 'C4')
                    pattern = 'Contracting 4-arm spiral';
                    patternnum = 16;
                elseif strcmpi(filename(1:2), 'C6')
                    pattern = 'Contracting 6-arm spiral';
                    patternnum = 17;
                else
                    warning('Unknown Pattern encountered');
                    keyboard;
                end
                
                setID = str2num(filename(3));
                
                split_filename = split(filename, '_');
                flightID = str2num(split_filename{2}(2:end));
                
                beeID = '99';
            else
                warning('Unknown directory encountered');
                keyboard;
            end 
        end
        
        
        
    end
end