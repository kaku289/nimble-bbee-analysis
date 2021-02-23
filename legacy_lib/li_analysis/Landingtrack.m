classdef Landingtrack < handle
    properties (Access = public)
        
        pattern = '';
        light = '';
        confirmedChecks_le = false; % If 2d and 3d checks are confirmed for this landing track

        nObjects = 0; % # of objects identified as landing tracks
        obj_id = []; 

        rawTrack = []; % 16 columns (obj_id	frame	timestamp	x	y	z	xvel	yvel	zvel P00	P01	P02	P11	P12	P22	P33	P44	P55)
        state = []; % Butterworth filtered state - 10 columns (timestamp	x	y	z	xvel	yvel	zvel xacc yacc zacc)
        
        scoredLandingParams = landingParameters.empty;
        
        % parameters
        time; % time vector for plots (t=0 at touchdown 1)
        displacement; % Eucledian distance from 3D position at touchdown 1
        
        state_le = [];
        state_t1 = [];
        state_t2 = [];
        
        duration_le_t1 = nan; % time_t1 - time_le
        duration_t1_t2 = nan; % time_t2 - time_t1
        duration_le_t2 = nan; % time_t2 - time_le
        
        dist_le_t1 = nan;
        disp_le_t1 = nan;
        turtuosity_le_t1 = nan;
        
        vbyz = nan;
        
        dt = 0.00571966171264648; % delta t between two consecutive frames (in s) from hdf5 files
        
        % parameters in Landing Disc Frame (LDF)
        state_LDF = [];
        state_le_LDF = [];
        state_t1_LDF = [];
        state_t2_LDF = [];
  
        equinox_state_LDF = []; % point where first oscillation is observed in landing trajectory
        stable_state_LDF = []; % traj between equinox_state and some y e.g. 0.15m - To get info about r_ref
        
        % For learning the parameters of the landing model
        track_subset_sysID = trackForLandingModel.empty;
                
    end
    methods
        function obj = Landingtrack(pattern, light, confirmedChecks, rawTrack, videoParams)
            if ~isempty(rawTrack)
                obj.pattern = pattern;
                obj.light = light;
                obj.confirmedChecks_le = confirmedChecks;

                obj.rawTrack = rawTrack;

                obj.obj_id = unique(rawTrack(:,1));
                obj.nObjects = length(obj.obj_id);

                obj.scoredLandingParams = copy(videoParams.landingParams);

                end_t = videoParams.timestamps(nanmax([obj.scoredLandingParams.legextension_framenumber ...
                        obj.scoredLandingParams.touchdown_framenumber ...
                        obj.scoredLandingParams.touchdown2_framenumber])) + 2*obj.dt;

                obj.rawTrack = sortrows(obj.rawTrack(obj.rawTrack(:,3)-(end_t) <= 1e-5,:), [3 2]);

                if ~isempty(obj.rawTrack) && sum(isnan(obj.rawTrack(:,3))) == 0
                    obj.filterTrack(20);
                end
            end
            
        end
        
        function filterTrack(obj, fc)
%             fc = 40;
            fs = 1/obj.dt; % fs is 174.8355 Hz
            % fs = 1/nanmean(diff(obj.rawTrack(:,3)));

            equidistant_time = [obj.rawTrack(1,3):obj.dt:obj.rawTrack(end,3)]';
            
            pos = interp1(obj.rawTrack(:,3)-obj.rawTrack(1,3), obj.rawTrack(:,4:6), equidistant_time-equidistant_time(1), 'makima');
            
            
            [num,den] = butter(2, fc/(fs/2),'low'); % 3rd order Butterworth filter %XXXX CHECK ORDER OF THE FILTER
            
            obj.state = zeros(size(pos,1),10);
            obj.state(:,1) = equidistant_time;
            obj.state(:,2:4) = filtfilt(num, den, pos);
            obj.state(:,5:7) =  diffxy(equidistant_time, obj.state(:,2:4)); %filtfilt(num, den, diffxy(equidistant_time, obj.state(:,2:4)));
            obj.state(:,8:10) = diffxy(equidistant_time, obj.state(:,2:4), [], 2); %filtfilt(num, den, diffxy(equidistant_time, obj.state(:,2:4), [], 2));
        end
        
        function compute_states_in_LDF(obj, landingDiscs)
            % Computes state in a inertial reference frame attached to the
            % center of the landing disc
            
            % landingDisc - array of instances of LandingDisc class
            landing_side = obj.scoredLandingParams.landing_side;
            landing_disc = landingDiscs(strcmpi({landingDiscs.side}, landing_side));
            
            assert(length(landing_disc) == 1, 'Multiple landing discs found! Not possible...');
            
            obj.state_LDF = Landingtrack.convert_to_landing_disc_reference_frame(obj.state, landing_disc.center', landing_side);
            obj.state_le_LDF = Landingtrack.convert_to_landing_disc_reference_frame(obj.state_le, landing_disc.center', landing_side);
            obj.state_t1_LDF = Landingtrack.convert_to_landing_disc_reference_frame(obj.state_t1, landing_disc.center', landing_side);
            obj.state_t2_LDF = Landingtrack.convert_to_landing_disc_reference_frame(obj.state_t2, landing_disc.center', landing_side);
            
        end
        
        function plotData(obj)
            
            % Plot trajectory
            trajPlot = figure;
            figure(trajPlot);
            figHandles(1) = trajPlot;
            subplotHandles(1) = subplot(3,1,1); hold on;
%             plot(obj.rawTrack(:,3)-obj.state_le(1),obj.rawTrack(:,4),'b.','MarkerSize',10);
            plot(obj.state(:,1)-obj.state_le(1),obj.state(:,2),'m.','MarkerSize',10);
            ylabel('x (m)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;
%             legend({'raw','Butterworth'});

            figHandles(end+1) = trajPlot;
            subplotHandles(end+1) = subplot(3,1,2); hold on;
%             plot(obj.rawTrack(:,3)-obj.state_le(1),obj.rawTrack(:,5),'b.','MarkerSize',10);
            plot(obj.state(:,1)-obj.state_le(1),obj.state(:,3),'m.','MarkerSize',10);
            ylabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;

            figHandles(end+1) = trajPlot;
            subplotHandles(end+1) = subplot(3,1,3); hold on;
%             plot(obj.rawTrack(:,3)-obj.state_le(1),obj.rawTrack(:,6),'b.','MarkerSize',10);
            plot(obj.state(:,1)-obj.state_le(1),obj.state(:,4),'m.','MarkerSize',10);
            ylabel('z (m)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;
            xlabel('time (s)', 'FontSize', 16);
            
            % Plot velocity
            velPlot = figure;
            figure(velPlot);
            figHandles(end+1) = velPlot;
            subplotHandles(end+1) = subplot(3,1,1); hold on;
%             plot(obj.rawTrack(:,3)-obj.state_le(1),obj.rawTrack(:,7),'b.','MarkerSize',10);
            plot(obj.state(:,1)-obj.state_le(1),obj.state(:,5),'m.','MarkerSize',10);
            ylabel('V_{gx} (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;
%             legend({'raw','Butterworth'});

            figHandles(end+1) = velPlot;
            subplotHandles(end+1) = subplot(3,1,2); hold on;
%             plot(obj.rawTrack(:,3)-obj.state_le(1),obj.rawTrack(:,8),'b.','MarkerSize',10);
            plot(obj.state(:,1)-obj.state_le(1),obj.state(:,6),'m.','MarkerSize',10);
            ylabel('V_{gy} (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;

            figHandles(end+1) = velPlot;
            subplotHandles(end+1) = subplot(3,1,3); hold on;
%             plot(obj.rawTrack(:,3)-obj.state_le(1),obj.rawTrack(:,9),'b.','MarkerSize',10);
            plot(obj.state(:,1)-obj.state_le(1),obj.state(:,7),'m.','MarkerSize',10);
            ylabel('V_{gz} (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;
            xlabel('time (s)', 'FontSize', 16);
            
            % Plot acceleration
            accPlot = figure;
            figure(accPlot);
            figHandles(end+1) = accPlot;
            subplotHandles(end+1) = subplot(3,1,1); hold on;
            plot(obj.state(:,1)-obj.state_le(1),obj.state(:,8),'m-o','MarkerSize',2,'LineWidth',1);
            ylabel('a_{gx} (m/s^2)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;
            legend({'Butterworth'});

            figHandles(end+1) = accPlot;
            subplotHandles(end+1) = subplot(3,1,2); hold on;
            plot(obj.state(:,1)-obj.state_le(1),obj.state(:,9),'m-o','MarkerSize',2,'LineWidth',1);
            ylabel('a_{gy} (m/s^2)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;

            figHandles(end+1) = accPlot;
            subplotHandles(end+1) = subplot(3,1,3); hold on;
            plot(obj.state(:,1)-obj.state_le(1),obj.state(:,10),'m-o','MarkerSize',2,'LineWidth',1);
            ylabel('a_{gz} (m/s^2)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;
            xlabel('time (s)', 'FontSize', 16);

            arrayfun(@(x,y) obj.drawlines(x,y), subplotHandles, figHandles);
            
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
            
            % Plot trajectory
            trajPlot = figure;
            figure(trajPlot);
            figHandles(1) = trajPlot;
            subplotHandles(1) = subplot(3,1,1); hold on;
            xaxis_dim{1} = 'time';
            plot(obj.state_LDF(:,1)-obj.state_le_LDF(1),obj.state_LDF(:,2),'m.','MarkerSize',10);
            ylabel('x (m)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;

            figHandles(end+1) = trajPlot;
            subplotHandles(end+1) = subplot(3,1,2); hold on;
            xaxis_dim{end+1} = 'time';
            plot(obj.state_LDF(:,1)-obj.state_le_LDF(1),obj.state_LDF(:,3),'m.','MarkerSize',10);
            ylabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;

            figHandles(end+1) = trajPlot;
            subplotHandles(end+1) = subplot(3,1,3); hold on;
            xaxis_dim{end+1} = 'time';
            plot(obj.state_LDF(:,1)-obj.state_le_LDF(1),obj.state_LDF(:,4),'m.','MarkerSize',10);
            ylabel('z (m)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;
            xlabel('time (s)', 'FontSize', 16);
            
            linkaxes(subplotHandles(1:3),'x');
            xlim([-1.5 inf]);
            
            % Plot velocity
            velPlot = figure;
            figure(velPlot);
            figHandles(end+1) = velPlot;
            subplotHandles(end+1) = subplot(3,1,1); hold on;
            xaxis_dim{end+1} = 'time';
            plot(obj.state_LDF(:,1)-obj.state_le_LDF(1),obj.state_LDF(:,5),'m.','MarkerSize',10);
            ylabel('V_{gx} (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;

            figHandles(end+1) = velPlot;
            subplotHandles(end+1) = subplot(3,1,2); hold on;
            xaxis_dim{end+1} = 'time';
            plot(obj.state_LDF(:,1)-obj.state_le_LDF(1),obj.state_LDF(:,6),'m.','MarkerSize',10);
            ylabel('V_{gy} (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;

            figHandles(end+1) = velPlot;
            subplotHandles(end+1) = subplot(3,1,3); hold on;
            xaxis_dim{end+1} = 'time';
            plot(obj.state_LDF(:,1)-obj.state_le_LDF(1),obj.state_LDF(:,7),'m.','MarkerSize',10);
            ylabel('V_{gz} (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;
            xlabel('time (s)', 'FontSize', 16);
            
            linkaxes(subplotHandles(4:6),'x');
            xlim([-1.5 inf]);
            
            % Plot acceleration
            accPlot = figure;
            figure(accPlot);
            figHandles(end+1) = accPlot;
            subplotHandles(end+1) = subplot(3,1,1); hold on;
            xaxis_dim{end+1} = 'time';
            plot(obj.state_LDF(:,1)-obj.state_le_LDF(1),obj.state_LDF(:,8),'m-o','MarkerSize',2,'LineWidth',1);
            ylabel('a_{gx} (m/s^2)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;

            figHandles(end+1) = accPlot;
            subplotHandles(end+1) = subplot(3,1,2); hold on;
            xaxis_dim{end+1} = 'time';
            plot(obj.state_LDF(:,1)-obj.state_le_LDF(1),obj.state_LDF(:,9),'m-o','MarkerSize',2,'LineWidth',1);
            ylabel('a_{gy} (m/s^2)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;

            figHandles(end+1) = accPlot;
            subplotHandles(end+1) = subplot(3,1,3); hold on;
            xaxis_dim{end+1} = 'time';
            plot(obj.state_LDF(:,1)-obj.state_le_LDF(1),obj.state_LDF(:,10),'m-o','MarkerSize',2,'LineWidth',1);
            ylabel('a_{gz} (m/s^2)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;
            xlabel('time (s)', 'FontSize', 16);
            
            linkaxes(subplotHandles(7:9),'x');
            xlim([-1.5 inf]);
            
            % Plot optical flow parameters
            % These parameters are plotted only for 1.5s seconds before LE
            state_subset = obj.state_LDF(obj.time>=-1.5,:);
            
            opticalExpansionPlot = figure;
            figure(opticalExpansionPlot);
            figHandles(end+1) = opticalExpansionPlot;
            subplotHandles(end+1) = subplot(2,1,1); hold on;
            xaxis_dim{end+1} = 'distance';
%             plot(obj.state_LDF(:,3),obj.state_LDF(:,6),'m.','MarkerSize',10);
            plot(state_subset(:,3),state_subset(:,6),'m.','MarkerSize',10);
            ylabel('V_{gy} (m/s)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;
            
            figHandles(end+1) = opticalExpansionPlot;
            subplotHandles(end+1) = subplot(2,1,2); hold on;
            xaxis_dim{end+1} = 'distance';
%             plot(obj.state_LDF(:,3),obj.state_LDF(:,6)./obj.state_LDF(:,3),'m.','MarkerSize',10);
            plot(state_subset(:,3),state_subset(:,6)./state_subset(:,3),'m.','MarkerSize',10);
            ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
            xlabel('y (m)', 'FontSize', 16);
            set(gca, 'FontSize', 18); grid on;
            
            linkaxes(subplotHandles(10:11),'x');
            
            arrayfun(@(x,y,z) obj.drawlines_LDF(x,y,z), subplotHandles, figHandles, xaxis_dim);
            
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
            ylabel('V_{gy} / y (1/s)', 'FontSize', 15);
%             ylabel('Optical rate of expansion (1/s)', 'FontSize', 15);
%             xlabel('Distance from the platform (m)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            ylim([-10 1]);
            
            figHandles(end+1) = opticalExpansionPlot;
            subplotHandles(end+1) = subplot(3,1,3); hold on;
            xaxis_dim{end+1} = 'distance';
%             plot(obj.state_LDF(:,3),obj.state_LDF(:,6)./obj.state_LDF(:,3),'m.','MarkerSize',10);
            plot(state_subset(:,3),state_subset(:,9)./state_subset(:,6),'b-','LineWidth',2,'MarkerSize',10);
            ylabel('a / V_{gy} (1/s)', 'FontSize', 15);
%             ylabel('Optical rate of expansion (1/s)', 'FontSize', 15);
%             xlabel('Distance from the platform (m)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            ylim([-50 50]);
            
% %             M = movvar(state_subset(:,6)./state_subset(:,3), [174 1]);
% %             figHandles(end+1) = opticalExpansionPlot;
% %             subplotHandles(end+1) = subplot(3,1,3); hold on;
% %             xaxis_dim{end+1} = 'distance';
% % %             plot(obj.state_LDF(:,3),obj.state_LDF(:,6)./obj.state_LDF(:,3),'m.','MarkerSize',10);
% %             plot(state_subset(:,3),M,'b.','MarkerSize',10);
% %             ylabel('Var. of r', 'FontSize', 15);
% %             xlabel('Distance from the platform (m)', 'FontSize', 15);
% %             set(gca, 'FontSize', 15); grid on;
% %             ylim([0 4]);
            
            linkaxes(subplotHandles(1:3),'x');
%             xlim([0 0.2]);
            
            arrayfun(@(x,y,z) obj.drawlines_Poster(x,y,z), subplotHandles, figHandles, xaxis_dim);
                        
            plotHandles = [opticalExpansionPlot];
        end
        
        function plotHandles = plotData_Time(obj)
            % Plot optical flow parameters
            % These parameters are plotted only for 1.5s seconds before LE
            
            str = '#375E98';
            color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

            state_subset = obj.state_LDF(obj.time>=-1.5,:);
            
            opticalExpansionPlot = figure;
            figure(opticalExpansionPlot);
            figHandles(1) = opticalExpansionPlot;
            subplotHandles(1) = subplot(5,1,1); hold on;
            xaxis_dim{1} = 'time';
            plot(state_subset(:,1)-obj.state_le_LDF(1,1),state_subset(:,6),'b.','MarkerSize',10);
            ylabel('V_{gy} (m/s)', 'FontSize', 15);
%             ylabel('Velocity towards the platform (m/s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            
            figHandles(end+1) = opticalExpansionPlot;
            subplotHandles(end+1) = subplot(5,1,2); hold on;
            xaxis_dim{end+1} = 'time';
            plot(state_subset(:,1)-obj.state_le_LDF(1,1),state_subset(:,3),'b.','MarkerSize',10);
            ylabel('y (m)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
%             ylim([-10 1]);
            
            figHandles(end+1) = opticalExpansionPlot;
            subplotHandles(end+1) = subplot(5,1,3); hold on;
            xaxis_dim{end+1} = 'time';
            plot(state_subset(:,1)-obj.state_le_LDF(1,1),state_subset(:,6)./state_subset(:,3),'b.','MarkerSize',10);
            ylabel('r (1/s)', 'FontSize', 15);
%             ylabel('Optical rate of expansion (1/s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
%             ylim([-10 1]);
            
            figHandles(end+1) = opticalExpansionPlot;
            subplotHandles(end+1) = subplot(5,1,4); hold on;
            xaxis_dim{end+1} = 'time';
            plot(state_subset(:,1)-obj.state_le_LDF(1,1),...
                state_subset(:,9)./state_subset(:,3)-state_subset(:,6).^2./state_subset(:,3).^2,'b.','MarkerSize',10);
            ylabel('r dot (1/s^2)', 'FontSize', 15);
            xlabel('Time to leg extension (s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
%             ylim([-10 1]);

            figHandles(end+1) = opticalExpansionPlot;
            subplotHandles(end+1) = subplot(5,1,5); hold on;
            xaxis_dim{end+1} = 'time';
            plot(state_subset(:,1)-obj.state_le_LDF(1,1),...
                state_subset(:,9),'b.','MarkerSize',10);
            ylabel('a_y (m/s^2)', 'FontSize', 15);
            xlabel('Time to leg extension (s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            
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
            xlim([-1.2 0]);
                        
            plotHandles = [opticalExpansionPlot];
        end
        
        function plotHandles = plotData_Time1(obj)
            % Plot optical flow parameters
            % These parameters are plotted only for 1.5s seconds before LE
            
            state_subset = obj.state_LDF(obj.time>=-1.5,:);
            
            opticalExpansionPlot = figure;
            figure(opticalExpansionPlot);
            figHandles(1) = opticalExpansionPlot;
           
            figHandles(end+1) = opticalExpansionPlot;
            subplotHandles(1) = subplot(3,1,1); hold on;
            xaxis_dim{1} = 'time';
            plot(state_subset(:,1)-obj.state_le_LDF(1,1),state_subset(:,6)./state_subset(:,3),'b','LineWidth',2,'MarkerSize',10);
            ylabel('V_{gy} / y (1/s)', 'FontSize', 15);
%             ylabel('Optical rate of expansion (1/s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
%             ylim([-10 1]);
            
            figHandles(end+1) = opticalExpansionPlot;
            subplotHandles(end+1) = subplot(3,1,2); hold on;
            xaxis_dim{end+1} = 'time';
            plot(state_subset(:,1)-obj.state_le_LDF(1,1),...
                state_subset(:,9)./state_subset(:,6),'b','LineWidth',2,'MarkerSize',10);
            ylabel('a_y / V_{gy} (1/s)', 'FontSize', 15);
            xlabel('Time to leg extension (s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            ylim([-80 80]);

            figHandles(end+1) = opticalExpansionPlot;
            subplotHandles(end+1) = subplot(3,1,3); hold on;
            xaxis_dim{end+1} = 'time';
            plot(state_subset(:,1)-obj.state_le_LDF(1,1),...
                state_subset(:,9),'b','LineWidth',2,'MarkerSize',10);
            ylabel('a_y (m/s^2)', 'FontSize', 15);
            xlabel('Time to leg extension (s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            
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
            xlim([-1.5 0]);
                        
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
        function state_LDF = convert_to_landing_disc_reference_frame(state, origin, landing_side)
            % state - N X 10 matrix similar to Landingtrack.state
            % origin - 1 X 3 vector containing origin of the landing
            % disc reference frame
            % landing_side - 'Hive' or 'Feeder'
            
            % returns state in landing disc reference frame (state_LDF)
                
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
    end
end