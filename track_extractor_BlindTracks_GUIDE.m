function varargout = track_extractor_BlindTracks_GUIDE(varargin)
% TRACK_EXTRACTOR_BLINDTRACKS_GUIDE MATLAB code for track_extractor_BlindTracks_GUIDE.fig
%      TRACK_EXTRACTOR_BLINDTRACKS_GUIDE, by itself, creates a new TRACK_EXTRACTOR_BLINDTRACKS_GUIDE or raises the existing
%      singleton*.
%
%      H = TRACK_EXTRACTOR_BLINDTRACKS_GUIDE returns the handle to a new TRACK_EXTRACTOR_BLINDTRACKS_GUIDE or the handle to
%      the existing singleton*.
%
%      TRACK_EXTRACTOR_BLINDTRACKS_GUIDE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACK_EXTRACTOR_BLINDTRACKS_GUIDE.M with the given input arguments.
%
%      TRACK_EXTRACTOR_BLINDTRACKS_GUIDE('Property','Value',...) creates a new TRACK_EXTRACTOR_BLINDTRACKS_GUIDE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before track_extractor_BlindTracks_GUIDE_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to track_extractor_BlindTracks_GUIDE_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help track_extractor_BlindTracks_GUIDE

% Last Modified by GUIDE v2.5 09-Jan-2020 11:11:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @track_extractor_BlindTracks_GUIDE_OpeningFcn, ...
                   'gui_OutputFcn',  @track_extractor_BlindTracks_GUIDE_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before track_extractor_BlindTracks_GUIDE is made visible.
function track_extractor_BlindTracks_GUIDE_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to track_extractor_BlindTracks_GUIDE (see VARARGIN)

% Choose default command line output for track_extractor_BlindTracks_GUIDE
handles.output = hObject;

% Loading tracks
try
    relevantTreatments = evalin("base","relevantTreatments");
    handles.tracks = [relevantTreatments.landingTracks]; % It is loading by reference.
    information = 'relevantTreatments variable exists! Loading it by reference...';
catch
    information = 'relevantTreatments variable does NOT exist :(';
end

% Loading appData
try
    handles.appData = evalin("base","appData"); % It is loading by reference.
    information = [information ' Also, loading existing appData.'];
catch
    handles.appData = appData_BlindtracksGUI();
    information = [information ' appData does NOT exist. Creating a new one.'];
end

% Loading info about relevant tracks
try
    ct_behaviour = evalin("base","ct_behaviour");
    ct_pattern = evalin("base","ct_pattern");
    ct_light = evalin("base","ct_light");
    behaviour = evalin("base","behaviour");
    pattern = evalin("base","pattern");
    light = evalin("base","light");
    
    handles.appData.filename = ['GUIDE_appData_Blindtracks_' pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour} '.mat'];
    information = [information ' Filename for saving appData: ' handles.appData.filename];
catch
    handles.appData.filename = ['GUIDE_appData_Blindtracks_' 'UNKNOWN' '.mat'];
    information = [information ' Filename for saving appData: ' handles.appData.filename];
end

%
handles.trackExcerpt = 1;

% Displaying info
displayInfo(hObject, handles, information);

% Setting axes properly
setAxes(hObject, handles);

% Setting mode
if handles.appData.mode == 1
    handles.mode2.Value = false;
    handles.mode1.Value = true;
elseif handles.appData.mode == 2
    handles.mode2.Value = true;
    handles.mode1.Value = false;
end

% Setting which part will be selected
handles.appData.part = 1;

% Setting states
% handles.gety_start = false;
% handles.gety_end = false;
% handles.y_start_xline = nan;
% handles.y_end_xline = nan;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes track_extractor_BlindTracks_GUIDE wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = track_extractor_BlindTracks_GUIDE_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function displayPartInfo(handles)
%  information - string to be displayed on displayBoard
handles.partDisplay.String = handles.appData.part;

function displayInfo(hObject, handles, information)
%  information - string to be displayed on displayBoard
handles.displayBoard.String = information;

function setAxes(hObject, handles)
% The first axes:
axes(handles.axes1); grid on;
ylabel('V_{gy} (m/s)', 'FontSize', 12);
% xlabel('y (m)', 'FontSize', 12);

% The second axes
axes(handles.axes2); grid on;
ylabel('V_{gy} / y (1/s)', 'FontSize', 12);
% xlabel('y (m)', 'FontSize', 12);

% The third axes
axes(handles.axes3); grid on;
ylabel('a_y / V_{gy} (1/s^2)', 'FontSize', 12);
xlabel('y (m)', 'FontSize', 12);

% The fourth axes
axes(handles.axes4); grid on; 
ylabel('V_{gy} / y (1/s)', 'FontSize', 12);
xlabel('y (m)', 'FontSize', 12);




function plotTrackData(hObject, handles, track)
% track - instance of BlindLandingtrack

part = handles.appData.part; % Part corresponds to which segment of mode 1or2 will be displayed

state_subset = track.state_LDF(handles.trackExcerpt).filteredState;

axes(handles.axes1);
plot(state_subset(:,3),state_subset(:,6),'m.','MarkerSize',10);
if handles.appData.mode == 2
    xline(handles.axes1,  track.DataGUI(handles.trackExcerpt).y(part,1), 'Label', 'y_{start}');
    xline(handles.axes1,  track.DataGUI(handles.trackExcerpt).y(part,2), 'Label', 'y_{end}');
elseif handles.appData.mode == 1
    xline(handles.axes1,  track.DataGUI(handles.trackExcerpt).y_mode1(part,1), 'Label', 'y_{start}');
    xline(handles.axes1,  track.DataGUI(handles.trackExcerpt).y_mode1(part,2), 'Label', 'y_{end}');
end
% grid on;

axes(handles.axes2);
plot(state_subset(:,3),state_subset(:,6)./state_subset(:,3),'m.','MarkerSize',10);
if handles.appData.mode == 2
    xline(handles.axes2,  track.DataGUI(handles.trackExcerpt).y(part,1), 'Label', 'y_{start}');
    xline(handles.axes2,  track.DataGUI(handles.trackExcerpt).y(part,2), 'Label', 'y_{end}');
elseif handles.appData.mode == 1
    xline(handles.axes2,  track.DataGUI(handles.trackExcerpt).y_mode1(part,1), 'Label', 'y_{start}');
    xline(handles.axes2,  track.DataGUI(handles.trackExcerpt).y_mode1(part,2), 'Label', 'y_{end}');
end
% grid on;

axes(handles.axes3);
% plot(state_subset(:,3),state_subset(:,9)./state_subset(:,3)-state_subset(:,6).^2./state_subset(:,3).^2,'m.','MarkerSize',10);
plot(state_subset(:,3),state_subset(:,9)./state_subset(:,6),'m','MarkerSize',10);

axes(handles.axes4);
hold off;
plot(track.state(handles.trackExcerpt).filteredState(:,3),track.state(handles.trackExcerpt).filteredState(:,6),'m.','MarkerSize',10); hold on;
plot(track.rawTrack(handles.trackExcerpt).rawState(:,5),track.rawTrack(handles.trackExcerpt).rawState(:,8),'b.','MarkerSize',10);

linkaxes([handles.axes1, handles.axes2, handles.axes3],'x');

xlim(handles.axes2, [min(state_subset(:,3))-0.01 0]);
ylim(handles.axes2, [-11 2]);
ylim(handles.axes3, [-40 40]);

setAxes(hObject, handles);

function updateGUIDisplay(hObject, handles)
part = handles.appData.part; % Part corresponds to which segment of mode 1or2 will be displayed
track = handles.tracks(handles.appData.currentTrack);

if handles.appData.mode == 2
    handles.y_start_editBox.String = num2str(track.DataGUI(handles.trackExcerpt).y(part,1), '%2.3f');
    handles.y_end_editBox.String = num2str(track.DataGUI(handles.trackExcerpt).y(part,2), '%2.3f');
elseif handles.appData.mode == 1
    handles.y_start_editBox.String = num2str(track.DataGUI(handles.trackExcerpt).y_mode1(part,1), '%2.3f');
    handles.y_end_editBox.String = num2str(track.DataGUI(handles.trackExcerpt).y_mode1(part,2), '%2.3f');
end

handles.checkbox1.Value = track.DataGUI(handles.trackExcerpt).isExcerptUseful;

% update part display
displayPartInfo(handles)





function displayBoard_Callback(hObject, eventdata, handles)
% hObject    handle to displayBoard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of displayBoard as text
%        str2double(get(hObject,'String')) returns contents of displayBoard as a double


% --- Executes during object creation, after setting all properties.
function displayBoard_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displayBoard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Loadstart.
function Loadstart_Callback(hObject, eventdata, handles)
% hObject    handle to Loadstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if handles.trackExcerpt == 0 %+ 1 <= length(track.state)
%     handles.trackExcerpt = handles.trackExcerpt + 1;
% end

information = ['Loaded track ' num2str(handles.appData.currentTrack) ' and excerpt ' num2str(handles.trackExcerpt)];
displayInfo(hObject, handles, information);
            
track = handles.tracks(handles.appData.currentTrack); % Instance of BlindLandingtrack
if isempty(track.DataGUI) || length(track.DataGUI) ~= length(track.state)
    for ct=1:length(track.state) % for each track excerpt
        track.DataGUI(ct) = DataGUI_BlindTracks();
    end
end


% To make backwards compatibility before 30/01/2020 - Comment it after
% processing high spokes first 800 trajectories
for ct=1:length(track.state) % for each track excerpt
    if size(track.DataGUI(ct).y_mode1,1) == 1 
        y_mode1 = track.DataGUI(ct).y_mode1;
        Vgy_mode1 = track.DataGUI(ct).Vgy_mode1;
        
        track.DataGUI(ct).y_mode1 = zeros(5,2);
        track.DataGUI(ct).y_mode1(1,:) = y_mode1;
        
        track.DataGUI(ct).Vgy_mode1 = zeros(5,2);
        track.DataGUI(ct).Vgy_mode1(1,:) = Vgy_mode1;
    end
end


plotTrackData(hObject, handles, track);
updateGUIDisplay(hObject, handles);

% Update handles structure
guidata(hObject, handles);


function y_start_editBox_Callback(hObject, eventdata, handles)
% hObject    handle to y_start_editBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_start_editBox as text
%        str2double(get(hObject,'String')) returns contents of y_start_editBox as a double


% --- Executes during object creation, after setting all properties.
function y_start_editBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_start_editBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_end_editBox_Callback(hObject, eventdata, handles)
% hObject    handle to y_end_editBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_end_editBox as text
%        str2double(get(hObject,'String')) returns contents of y_end_editBox as a double


% --- Executes during object creation, after setting all properties.
function y_end_editBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_end_editBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, event, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
key = event.Key;
switch handles.zoomin.State
    case 'off'
        if strcmpi(key, 'numpad1')
            part = handles.appData.part;
            
            displayInfo(hObject, handles, 'In getting y_start mode... Double click to get a point in first axes.');
            
            [xi, yi] = getpts(handles.axes1);
            track = handles.tracks(handles.appData.currentTrack);
            
            % Snap onto the nearest point in the dataset
            % Assumes that point is acquired in axes1 (V vs y)
            state_subset = track.state_LDF(handles.trackExcerpt).filteredState;
            [~,indx] = min(sqrt((state_subset(:,3)-xi(end)).^2 +(state_subset(:,6)-yi(end)).^2));
            if handles.appData.mode == 2
                track.DataGUI(handles.trackExcerpt).y(part,1) = state_subset(indx,3);
                track.DataGUI(handles.trackExcerpt).Vgy(part,1) = state_subset(indx,6);
            elseif handles.appData.mode == 1
                track.DataGUI(handles.trackExcerpt).y_mode1(part,1) = state_subset(indx,3);
                track.DataGUI(handles.trackExcerpt).Vgy_mode1(part,1) = state_subset(indx,6);
            end
            handles.y_start_editBox.String = num2str(state_subset(indx,3), '%2.3f');
            
            for ct=1:length(handles.axes2.Children)
                if isa(handles.axes2.Children(ct),'matlab.graphics.chart.decoration.ConstantLine') && ...
                        strcmpi(handles.axes2.Children(ct).Label, 'y_{start}')
                    delete(handles.axes2.Children(ct));
                    break;
                end
            end
            xline(handles.axes2, state_subset(indx,3), 'Label', 'y_{start}');
            
            for ct=1:length(handles.axes1.Children)
                if isa(handles.axes1.Children(ct),'matlab.graphics.chart.decoration.ConstantLine') && ...
                        strcmpi(handles.axes1.Children(ct).Label, 'y_{start}')
                    delete(handles.axes1.Children(ct));
                    break;
                end
            end
            xline(handles.axes1, state_subset(indx,3), 'Label', 'y_{start}');
            
%             if handles.appData.mode == 2
%                 track.DataGUI(handles.trackExcerpt).y(1,1) = xi(end-1);
%                 track.DataGUI(handles.trackExcerpt).Vgy(1,1) = yi(end-1);
%             elseif handles.appData.mode == 1
%                 track.DataGUI(handles.trackExcerpt).y_mode1(1,1) = xi(end-1);
%                 track.DataGUI(handles.trackExcerpt).Vgy_mode1(1,1) = yi(end-1);
%             end
%             
%             handles.y_start_editBox.String = num2str(xi(end-1), '%2.3f');
%             
%             for ct=1:length(handles.axes2.Children)
%                 if isa(handles.axes2.Children(ct),'matlab.graphics.chart.decoration.ConstantLine') && ...
%                         strcmpi(handles.axes2.Children(ct).Label, 'y_{start}')
%                     delete(handles.axes2.Children(ct));
%                     break;
%                 end
%             end
%             xline(handles.axes2, xi(end-1), 'Label', 'y_{start}');
%             
%             for ct=1:length(handles.axes1.Children)
%                 if isa(handles.axes1.Children(ct),'matlab.graphics.chart.decoration.ConstantLine') && ...
%                         strcmpi(handles.axes1.Children(ct).Label, 'y_{start}')
%                     delete(handles.axes1.Children(ct));
%                     break;
%                 end
%             end
%             xline(handles.axes1, xi(end-1), 'Label', 'y_{start}');
            
            displayInfo(hObject, handles, '');
            
            % Update handles structure
            guidata(hObject, handles);
            
        elseif strcmpi(key, 'numpad2')
            part = handles.appData.part;
            
            displayInfo(hObject, handles, 'In getting y_end mode... Double click to get a point in bottom figure.');
            
            [xi, yi] = getpts(handles.axes1);
            track = handles.tracks(handles.appData.currentTrack);
            % Snap onto the nearest point in the dataset
            % Assumes that point is acquired in axes1 (V vs y)
            state_subset = track.state_LDF(handles.trackExcerpt).filteredState;
            [~,indx] = min(sqrt((state_subset(:,3)-xi(end)).^2 +(state_subset(:,6)-yi(end)).^2));
            if handles.appData.mode == 2
                track.DataGUI(handles.trackExcerpt).y(part,2) = state_subset(indx,3);
                track.DataGUI(handles.trackExcerpt).Vgy(part,2) = state_subset(indx,6);
            elseif handles.appData.mode == 1
                track.DataGUI(handles.trackExcerpt).y_mode1(part,2) = state_subset(indx,3);
                track.DataGUI(handles.trackExcerpt).Vgy_mode1(part,2) = state_subset(indx,6);
            end
            
            handles.y_end_editBox.String = num2str(state_subset(indx,3), '%2.3f');
            
            for ct=1:length(handles.axes2.Children)
                if isa(handles.axes2.Children(ct),'matlab.graphics.chart.decoration.ConstantLine') && ...
                        strcmpi(handles.axes2.Children(ct).Label, 'y_{end}')
                    delete(handles.axes2.Children(ct));
                    break;
                end
            end
            xline(handles.axes2, state_subset(indx,3), 'Label', 'y_{end}');
            
            for ct=1:length(handles.axes1.Children)
                if isa(handles.axes1.Children(ct),'matlab.graphics.chart.decoration.ConstantLine') && ...
                        strcmpi(handles.axes1.Children(ct).Label, 'y_{end}')
                    delete(handles.axes1.Children(ct));
                    break;
                end
            end
            xline(handles.axes1, state_subset(indx,3), 'Label', 'y_{end}');            
            
            
%             if handles.appData.mode == 2
%                 track.DataGUI(handles.trackExcerpt).y(1,2) = xi(end-1);
%                 track.DataGUI(handles.trackExcerpt).Vgy(1,2) = yi(end-1);
%             elseif handles.appData.mode == 1
%                 track.DataGUI(handles.trackExcerpt).y_mode1(1,2) = xi(end-1);
%                 track.DataGUI(handles.trackExcerpt).Vgy_mode1(1,2) = yi(end-1);
%             end
%             
%             handles.y_end_editBox.String = num2str(xi(end-1), '%2.3f');
%             
%             for ct=1:length(handles.axes2.Children)
%                 if isa(handles.axes2.Children(ct),'matlab.graphics.chart.decoration.ConstantLine') && ...
%                         strcmpi(handles.axes2.Children(ct).Label, 'y_{end}')
%                     delete(handles.axes2.Children(ct));
%                     break;
%                 end
%             end
%             xline(handles.axes2, xi(end-1), 'Label', 'y_{end}');
%             
%             for ct=1:length(handles.axes1.Children)
%                 if isa(handles.axes1.Children(ct),'matlab.graphics.chart.decoration.ConstantLine') && ...
%                         strcmpi(handles.axes1.Children(ct).Label, 'y_{end}')
%                     delete(handles.axes1.Children(ct));
%                     break;
%                 end
%             end
%             xline(handles.axes1, xi(end-1), 'Label', 'y_{end}');
            
            displayInfo(hObject, handles, '');
            
            % Update handles structure
            guidata(hObject, handles);
            
        elseif strcmpi(key, 'numpad4')
            switch handles.zoomin.State
                case 'off'
                    zoom(handles.axes2, 'on');
                    
                    % Re-enable keypress capture in pan or zoom mode
                    % http://undocumentedmatlab.com/blog/enabling-user-callbacks-during-zoom-pan
                    hManager = uigetmodemanager(handles.figure1);
                    [hManager.WindowListenerHandles.Enabled] = deal(false);
                    
                    set(handles.figure1, 'KeyPressFcn', @(hObject,eventdata)track_extractor_BlindTracks_GUIDE('figure1_KeyPressFcn',hObject,event,guidata(hObject)));
                    %             set(hFig, 'KeyPressFcn', @figure1_KeyPressFcn);
                case 'on'
                    zoom(handles.axes2, 'off');
            end
        elseif strcmpi(key, 'numpad5')
            zoom(handles.axes2,'out');
        elseif strcmpi(key, 'leftarrow')
            prevTrackButton_Callback(hObject, event, handles);
        elseif strcmpi(key, 'rightarrow')
            nextTrackButton_Callback(hObject, event, handles);
        elseif strcmpi(key, 'uparrow')
            saveDataButton_Callback(hObject, event, handles);
        elseif strcmpi(key, 'numpad7') % toggle checkbox1 (whether track is useful for sysID or not)
            track = handles.tracks(handles.appData.currentTrack);
            track.DataGUI(handles.trackExcerpt).isExcerptUseful = ~track.DataGUI(handles.trackExcerpt).isExcerptUseful;
            handles.checkbox1.Value = track.DataGUI(handles.trackExcerpt).isExcerptUseful;
            
            % Update handles structure
            guidata(hObject, handles);
        elseif strcmpi(key, 'numpad8') % toggle mode of the app
            if handles.appData.mode == 1
                mode2_Callback(hObject, nan, handles);
            elseif handles.appData.mode == 2
                mode1_Callback(hObject, nan, handles)
            end
        elseif strcmpi(key, 'numpad3')
            % Switching between "part"
            if rem(handles.appData.part+1,3) ~= 0
                handles.appData.part = rem(handles.appData.part+1,3);
            else
                handles.appData.part = 3;
            end
            
            % Update handles structure
            guidata(hObject, handles);
            
            % Load new data
            Loadstart_Callback(hObject, event, handles);          
            
        end
        
    case 'on'
        zoom(handles.axes2, 'off');
end

% if handles.gety_start
%     
% end




% --- Executes on button press in resetButton.
function resetButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if handles.trackExcerpt == 0
%     handles.trackExcerpt = 1;
% end

part = handles.appData.part; % Part corresponds to part-th segment of mode 1or2 
track = handles.tracks(handles.appData.currentTrack);
if handles.appData.mode == 2
    track.DataGUI(handles.trackExcerpt).y(part,1) = 0;
    track.DataGUI(handles.trackExcerpt).y(part,2) = 0;
elseif handles.appData.mode == 1
    track.DataGUI(handles.trackExcerpt).y_mode1(part,1) = 0;
    track.DataGUI(handles.trackExcerpt).y_mode1(part,2) = 0;
end

indx2del = [];
for ct=1:length(handles.axes2.Children)
    if isa(handles.axes2.Children(ct),'matlab.graphics.chart.decoration.ConstantLine')
        if strcmpi(handles.axes2.Children(ct).Label, 'y_{end}') || strcmpi(handles.axes2.Children(ct).Label, 'y_{start}')
            indx2del = [indx2del; ct];
        end
    end
end
delete(handles.axes2.Children(indx2del));
xline(handles.axes2, 0, 'Label', 'y_{start}');
xline(handles.axes2, 0, 'Label', 'y_{end}');

indx2del = [];
for ct=1:length(handles.axes1.Children)
    if isa(handles.axes1.Children(ct),'matlab.graphics.chart.decoration.ConstantLine')
        if strcmpi(handles.axes1.Children(ct).Label, 'y_{end}') || strcmpi(handles.axes1.Children(ct).Label, 'y_{start}')
            indx2del = [indx2del; ct];
        end
    end
end
delete(handles.axes1.Children(indx2del));
xline(handles.axes1, 0, 'Label', 'y_{start}');
xline(handles.axes1, 0, 'Label', 'y_{end}');

updateGUIDisplay(hObject, handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in nextTrackButton.
function nextTrackButton_Callback(hObject, eventdata, handles)
% hObject    handle to nextTrackButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
track = handles.tracks(handles.appData.currentTrack);
if handles.appData.currentTrack < length(handles.tracks)
    
    if handles.trackExcerpt == 0 || handles.trackExcerpt < length(track.state_LDF)
        handles.trackExcerpt = handles.trackExcerpt + 1;
    else
        handles.trackExcerpt = 1;
        handles.appData.currentTrack = handles.appData.currentTrack + 1;
    end
    
    % Setting the part
    handles.appData.part = 1; % Part corresponds to which segment of mode 1or2 will be displayed
    
    % Update handles structure
    guidata(hObject, handles);
    
    Loadstart_Callback(hObject, eventdata, handles);
%     keyboard
%     track = handles.tracks(handles.appData.currentTrack);
%     plotTrackData(hObject, handles, track);
%     updateGUIDisplay(hObject, handles)
else
    displayInfo(hObject, handles, 'Can not go beyond this track!');
end


% --- Executes on button press in prevTrackButton.
function prevTrackButton_Callback(hObject, eventdata, handles)
% hObject    handle to prevTrackButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Setting the part
handles.appData.part = 1; % Part corresponds to which segment of mode 1or2 will be displayed
    
if handles.appData.currentTrack==1 && handles.trackExcerpt==1
    displayInfo(hObject, handles, 'Can not go beyond this track!');

elseif handles.appData.currentTrack==1 && handles.trackExcerpt > 1
    handles.trackExcerpt = handles.trackExcerpt - 1;
    
    % Update handles structure
    guidata(hObject, handles);
    
    Loadstart_Callback(hObject, eventdata, handles);

elseif handles.appData.currentTrack > 1
    
    if handles.trackExcerpt > 1
        handles.trackExcerpt = handles.trackExcerpt - 1;
    else
        handles.appData.currentTrack = handles.appData.currentTrack - 1;
        track = handles.tracks(handles.appData.currentTrack);
        handles.trackExcerpt = length(track.state_LDF);
    end
    
    % Update handles structure
    guidata(hObject, handles);
    
    Loadstart_Callback(hObject, eventdata, handles);
    
%     track = handles.tracks(handles.appData.currentTrack);
%     plotTrackData(hObject, handles, track);
%     updateGUIDisplay(hObject, handles)
end


% --- Executes on button press in saveDataButton.
function saveDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
displayInfo(hObject, handles, 'Saving data! Wait for 2 minutes..');

% Update handles structure
guidata(hObject, handles);

pause(0.1);

treatments = evalin("base","treatments");
appData = handles.appData;

save(fullfile(handles.appData.dataDir, handles.appData.filename),'appData');
save(fullfile(handles.appData.dataDir,'BlindLandingtracks.mat'),'treatments');

displayInfo(hObject, handles, 'Data Saved!');
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in mode1.
function mode1_Callback(hObject, eventdata, handles)
% hObject    handle to mode1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mode1
handles.appData.mode = 1;
handles.mode2.Value = false;
handles.mode1.Value = true;
% Update handles structure
guidata(hObject, handles);
Loadstart_Callback(hObject, eventdata, handles)

% --- Executes on button press in mode2.
function mode2_Callback(hObject, eventdata, handles)
% hObject    handle to mode2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mode2
handles.appData.mode = 2;
handles.mode2.Value = true;
handles.mode1.Value = false;
% Update handles structure
guidata(hObject, handles);
Loadstart_Callback(hObject, eventdata, handles);
% keyboard;



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
% keyboard;
handles.appData.currentTrack = str2num(handles.edit4.String{:});
handles.trackExcerpt = 1;
% Update handles structure
guidata(hObject, handles);
Loadstart_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
