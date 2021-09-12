clc; close all;
% clear;
% 
% data from data_write used in statistical analysis of rdot

%% Extract rrefEntry data for statistical analysis
clc; close all;
treatments = treatments(1:14*8); % Taking experiments for first 14 days
% combination

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high', 'lower'};
behaviour = {'rising','constant','sleeping'};

% data = data4rrefEntry.empty;
hasTakeoff = []; % whether or not the landing track has takeoff dynamics in it as defined in filteredState_BlindLandingtrack.hasTakeoff(..) method

rref_median = 2.87; rref_bin = 1;
y0_median = 0.21; y0_bin = 0.06;
chosen_fac = 1.5;

data = struct.empty;
clear dummy;
figure;
for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        for ct_behaviour = 2%1:length(behaviour)
        
           % Selecting relevant treatments
            if strcmpi(behaviour{ct_behaviour}, 'rising')
                relevantTreatments = treatments(strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
                                     strcmpi({treatments.light}, light{ct_light}) & ...
                                     rem(1:length(treatments), 8)==1);
            elseif strcmpi(behaviour{ct_behaviour}, 'constant')
                relevantTreatments = treatments(strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
                                     strcmpi({treatments.light}, light{ct_light}) & ...
                                     rem(1:length(treatments), 8)>1 & ...
                                     rem(1:length(treatments), 8)<8);
            elseif strcmpi(behaviour{ct_behaviour}, 'sleeping')
                relevantTreatments = treatments(strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
                                     strcmpi({treatments.light}, light{ct_light}) & ...
                                     rem(1:length(treatments), 8)==0);
            else
                error('What other treatments did you perform dude?')
            end
            
            
            for ct_treatment=1:length(relevantTreatments) % for each relevant treatment
                treatment = relevantTreatments(ct_treatment);
                stateLDF = [treatment.landingTracks.state_LDF]';
                
                arrayfun(@(x) x.setLandingSide(),[treatment.landingTracks.state_LDF]'); % to store landing side in the rrefSegments
                
                data1 = arrayfun(@(x) x.rrefEntrySegments,[treatment.landingTracks.state_LDF]','UniformOutput',false);
                data1 = horzcat(data1{:});
                
                
                
                hastakeoff_pertreatment = arrayfun(@(x) x.hasTakeoff(treatment.landingDiscs)*ones(1,length(x.rrefEntrySegments)),[treatment.landingTracks.state_LDF]','UniformOutput',false);
                hastakeoff_pertreatment = horzcat(hastakeoff_pertreatment{:});
                
                stateLDFindx = arrayfun(@(x) x*ones(1,length(stateLDF(x).rrefEntrySegments)),1:length(stateLDF),'UniformOutput',false);
                stateLDFindx = horzcat(stateLDFindx{:});
                
                % Discard empty intervals
%                 indices = arrayfun(@(x) ~isempty(x.intervals) & any(-x.rref>=rref_median-rref_bin/2) & any(-x.rref<=rref_median+rref_bin/2) & ...
%                                     any(-x.yEntryStart>=y0_median-y0_median/2) & any(-x.yEntryStart<=y0_median+y0_median/2), data1);      
                indices = arrayfun(@(x) ~isempty(x.intervals) & abs(x.factor-chosen_fac)<1e-6, data1); 
                % 
%                 if any(indices)
%                     keyboard;
%                 end
                
                % useful data
                rrefEntrySegments = data1(indices); % non-empty rrefEntrySegments
                stateLDFindx = stateLDFindx(indices);
                hastakeoff = hastakeoff_pertreatment(indices);
                
                for ct=1:length(rrefEntrySegments)
%                     if hastakeoff(ct)
                        entrySegment = rrefEntrySegments(ct);
                        track = stateLDF(stateLDFindx(ct));
                        
                        x = entrySegment;
%                         indices = (-x.rref>=rref_median-rref_bin/2) & (-x.rref<=rref_median+rref_bin/2) & ...
%                                     (-x.yEntryStart>=y0_median-y0_median/2) & (-x.yEntryStart<=y0_median+y0_median/2);
                        indices = true(length(x.rref),1);
                        
                        intervals = x.intervals(indices,:);
                        rref = x.rref(indices);
                        acc_actual = x.acc_actual(indices);
                        acc_rdotsim = x.acc_rdotsim(indices);
                        
                        
                        
                        for ct1=1:length(rref)
%                             dummy.a = track.filteredState(intervals(ct1,1):intervals(ct1,2),9); 
%                             a = acc_actual{ct1}; % actual accleration during entry
                            dummy.a = acc_rdotsim{ct1}; % acc from rdot sim during entry
                            dummy.deltar = track.filteredState(intervals(ct1,1):intervals(ct1,2),6)./track.filteredState(intervals(ct1,1):intervals(ct1,2),3) - rref(ct1);
                            dummy.rref = -rref(ct1);
                            dummy.track = track;
                            dummy.y0 = x.yEntryStart(ct1);
                            dummy.delta_r0 = x.delta_r(ct1);
                            dummy.hasTakeoff = hastakeoff(ct);
                            dummy.light = ct_light;
                            dummy.amean = x.delta_Ventry(ct1)/x.delta_tentry(ct1);
                            [dummy.R,dummy.P] = corrcoef(dummy.deltar, dummy.a);
                            data = [data; dummy];
                            
%                             plot(dummy.deltar, dummy.a, 'b'); hold on;
%                             plot(dummy.deltar, a, 'r')
%                             pause(0.5)
%                             close all
                            
                        end
%                     end
                end               
               
%                 size(data)
%                 keyboard;
            end
        end
    end
end

%%
a = vertcat(data.a); deltar = vertcat(data.deltar); 
% figure;
% plot(deltar, a, '.')

figure;
corr_val = arrayfun(@(x) x.R(1,2), data);
p_val = arrayfun(@(x) x.P(1,2), data);
subplot(1,2,1);
histogram(corr_val);
subplot(1,2,2);
histogram(p_val);

[sum(p_val<0.05) sum(corr_val<0) sum(p_val<0.05 & corr_val<0) length(p_val)]./length(p_val)
% createBoxPlot({corr_val, p_val}, {'Corr coef', 'p-value'}', 'Correlation between A(t) and deltar(t)')

figure;
histogram(a); 

figure;
histogram(deltar); 



%% light=3, from free-flight [sum(p_val<0.05) sum(corr_val<0) sum(p_val<0.05 & corr_val<0) length(p_val)]./length(p_val)
% 0.5677    0.7398    0.4657    1.0000

% Whole data
% 0.6077    0.7625    0.5157    1.0000

%% Plot supplementary figure
clc; close all;
map = brewermap(3,'Set1'); 

data1 = data([data.amean]>0);

a = vertcat(data1.a);

hastakeoff = arrayfun(@(x) x.hasTakeoff*ones(length(x.a),1),data1,'UniformOutput',false);
hastakeoff = vertcat(hastakeoff{:});

light = arrayfun(@(x) x.light*ones(length(x.a),1),data1,'UniformOutput',false);
light = vertcat(light{:});

rref = arrayfun(@(x) x.rref*ones(length(x.a),1),data1,'UniformOutput',false);
rref = vertcat(rref{:});

delta_r0 = arrayfun(@(x) x.delta_r0*ones(length(x.a),1),data1,'UniformOutput',false);
delta_r0 = vertcat(delta_r0{:});

y0 = arrayfun(@(x) x.y0*ones(length(x.a),1),data1,'UniformOutput',false);
y0 = vertcat(y0{:});

figure;
histogram(a(hastakeoff==1),'facecolor',map(2,:),'facealpha',.5,'edgecolor','none','Normalization','pdf');
hold on;
histogram(a(hastakeoff==0),'facecolor',map(3,:),'facealpha',.5,'edgecolor','none','Normalization','pdf');
legend('From take-off','From free-flight','fontsize',16)
xlabel('Instantaneous acc (ms-2)', 'FontSize', 16);
ylabel('Probability density', 'FontSize', 16);
set(gca, 'FontSize', 16);
xlim([-3 6]);

figure;
histogram(a(light==1),'facecolor',map(2,:),'facealpha',.5,'edgecolor','none','Normalization','pdf');
hold on;
histogram(a(light==3),'facecolor',map(3,:),'facealpha',.5,'edgecolor','none','Normalization','pdf');
legend('low light','high light','fontsize',16)
xlabel('Instantaneous acc (ms-2)', 'FontSize', 16);
ylabel('Probability density', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure;
histogram2(rref,a,'DisplayStyle','tile','Normalization','pdf');
ylabel('Instantaneous acc (ms-2)', 'FontSize', 16);
xlabel('rref (s-1)', 'FontSize', 16);
set(gca, 'FontSize', 16);
colorbar

figure;
histogram2(delta_r0,a,'DisplayStyle','tile','Normalization','pdf');
ylabel('Instantaneous acc (ms-2)', 'FontSize', 16);
xlabel('delta r (s-1)', 'FontSize', 16);
set(gca, 'FontSize', 16);
colorbar

figure;
histogram2(y0,a,'DisplayStyle','tile','Normalization','pdf');
ylabel('Instantaneous acc (ms-2)', 'FontSize', 16);
xlabel('y0 (m)', 'FontSize', 16);
set(gca, 'FontSize', 16);
colorbar
% xlim([-3 6]);

%% %
