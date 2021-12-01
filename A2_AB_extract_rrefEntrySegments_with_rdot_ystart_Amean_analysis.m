%% Plot instantaneous acclerations with different variables (supplement figure)
% load data_write for r file
clc;
close all;

chosen_fac = 1;
approachCol = 1;
landingSideCol = 2;
windCol = 3;
timeCol = 4;
dayCol = 5;
yCol = 13;
rrefCol = 7;
rdotCol = 8;
factorCol = 9;
isriseCol = 10;
deltarCol = 11;
deltaVCol = 16;
deltatCol = 17;
hasTakeoffCol = 12;
ameanCol = 18;
data_write(:,18) = data_write(:,deltaVCol)./data_write(:,deltatCol);


% factors = [0.25:0.25:2.5];
factors = [1];
for ct_factor=1:length(factors)
    factor = factors(ct_factor);
    
    data_fac = data(abs([data.factor]-factor)<1e-6)';
    hasTakeoff_fac = hasTakeoff(abs([data.factor]-factor)<1e-6)';
    N = length(data_fac);
    
    % Create nominal and ordinal variables
%     approach = arrayfun(@(i) i*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
%     side = arrayfun(@(i) data_fac(i).side*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
%     pattern = arrayfun(@(i) data_fac(i).pattern*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false) ;
%     light = arrayfun(@(i) data_fac(i).light*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
%     time = arrayfun(@(i) data_fac(i).time*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
%     day = arrayfun(@(i) data_fac(i).day*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
    
%     hasTakeoff_fac = arrayfun(@(i) hasTakeoff_fac(i)*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
    
    
end
% Plotting acclerations for hastakeoff factor
acc_actual_hastakeoff0 = vertcat(data_fac(hasTakeoff_fac==0).acc_actual);
acc_actual_hastakeoff0 = vertcat(acc_actual_hastakeoff0{:});
acc_actual_hastakeoff1 = vertcat(data_fac(hasTakeoff_fac==1).acc_actual);
acc_actual_hastakeoff1 = vertcat(acc_actual_hastakeoff1{:});
map = brewermap(3,'Set1'); 
figure;
histogram(acc_actual_hastakeoff1,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none','Normalization','pdf');
hold on;
histogram(acc_actual_hastakeoff0,'facecolor',map(3,:),'facealpha',.5,'edgecolor','none','Normalization','pdf');
legend('From take-off','From free-flight','fontsize',16)
xlabel('acc actual (m/s-2)', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);
xlim([-4 8]);

acc_actual_light1 = vertcat(data_fac([data_fac.light] == 1).acc_actual);
acc_actual_light1 = vertcat(acc_actual_light1{:});
acc_actual_light3 = vertcat(data_fac([data_fac.light] == 3).acc_actual);
acc_actual_light3 = vertcat(acc_actual_light3{:});
map = brewermap(3,'Set1'); 
figure;
histogram(acc_actual_light1,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none','Normalization','pdf');
hold on;
histogram(acc_actual_light3,'facecolor',map(3,:),'facealpha',.5,'edgecolor','none','Normalization','pdf');
legend('low light','high light','fontsize',16)
xlabel('acc actual (m/s-2)', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);
xlim([-4 8]);


% acc_rdotsim_hastakeoff0 = vertcat(data_fac(hasTakeoff_fac==0).acc_rdotsim);
% acc_rdotsim_hastakeoff0 = vertcat(acc_rdotsim_hastakeoff0{:});
% acc_rdotsim_hastakeoff1 = vertcat(data_fac(hasTakeoff_fac==1).acc_rdotsim);
% acc_rdotsim_hastakeoff1 = vertcat(acc_rdotsim_hastakeoff1{:});
% map = brewermap(3,'Set1'); 
% figure;
% histogram(acc_rdotsim_hastakeoff1,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none');
% hold on;
% histogram(acc_rdotsim_hastakeoff0,'facecolor',map(3,:),'facealpha',.5,'edgecolor','none');
% legend('From take-off','From free-flight','fontsize',16)
% xlabel('acc rdotsim (m/s-2)', 'FontSize', 16);
% ylabel('Occurences', 'FontSize', 16);
% set(gca, 'FontSize', 16);
% axis tight

%% Plot histogram of Amean
chosen_fac = 1;
% dayCol = 6;
% rdotCol =9;
% factorCol = 10;
% isriseCol = 11;
ameanCol = 18;
% data_ss = data_write(data_write(:,isriseCol) == 1 & data_write(:,ameanCol) >0 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6, :);
data_ss = data_write(data_write(:,isriseCol) == 1 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6, :);

figure;
histogram(data_ss(:,ameanCol));
% xlim([-Inf 5.5]);
% xticks([-1:1:5]);
% histogram(-vertcat(data.rmean), [0:0.5:8]);
% histogram(-vertcat(data.rmean), [0:0.5:9.5]);

xlabel('Mean accleration during entry segments, Amean (ms-2)', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);
% pd = fitdist(data_ss(:,rdotCol),'Gamma');
median(data_ss(:,ameanCol))

% data used in analysis
data_ss = data_write(data_write(:,isriseCol) == 1 & data_write(:,ameanCol) >0 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6, :);

%% Plot results from Amean model
%% Plot statistical model from R
close all;
% R model

% Factor = 1, with ystart as y
coeffs = [1.97341442  0.52876583  0.17191475 0.39327473 0.59033427 -1.63244769 -0.67882068 -0.10296793];
% coeffs = [1.58702290  0.42800399  0.05374526  0.16907431  0.32467186 -1.45318255 -0.58756146 0];

% Factor = 1.5, with ystart as y and inclduing hasTakeoff as a factor
% coeffs = [1.77080369  0.50743000  0.04233314  0.14084481  0.42358406 -1.47303718 -0.55294707 0]; % [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]

% data_all.approach = data_write(:,1);
% data_all.landingSide = data_write(:,2);
% data_all.pattern = data_write(:,3);
% data_all.light = data_write(:,4);
% data_all.day = data_write(:,6);
% data_all.y = data_write(:,7);
% data_all.rref = data_write(:,8);
% data_all.rdot = data_write(:,9);
% data_all.factor = data_write(:,10);
% data_all.isRise = data_write(:,11);
approachCol = 1;
landingSideCol = 2;
windCol = 3;
timeCol = 4;
dayCol = 5;
yCol = 13;
rrefCol = 7;
rdotCol = 8;
factorCol = 9;
isriseCol = 10;
deltarCol = 11;
deltaVCol = 16;
deltatCol = 17;
hasTakeoffCol = 12;
% data_write(:,end+1) = data_write(:,rrefCol)-data_write(:,deltarCol);
% r0Col = size(data_write,2);
data_write(:,18) = data_write(:,deltaVCol)./data_write(:,deltatCol);
dvdtCol = size(data_write,2);
% rdotCol = ameanCol; %dvdtCol;
ameanCol = 18;

% get subset of data

chosen_fac = 1;
% data_ss = data_write(abs(data_write(:,factorCol) - chosen_fac) < 1e-6 & data_write(:,dayCol) <= 20190710, :);
% data_ss = data_write(data_write(:,isriseCol) == 1 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6 & data_write(:,dayCol) <= 20190710, :);
data_ss = data_write(data_write(:,isriseCol) == 1 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6 & data_write(:,dayCol) <= 20190710 & data_write(:,ameanCol) > 0, :);
data_ss = data_write(data_write(:,isriseCol) == 1 & data_write(:,hasTakeoffCol) == 0 & data_write(:,ameanCol) >0 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6, :);
% data_ss = data_write(data_write(:,isriseCol) == 1 & data_write(:,hasTakeoffCol) == 1 & data_write(:,ameanCol) >0 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6 & data_write(:,dayCol) <= 20190710, :);

% Create basic plots
figure;
subplot(1,2,1)
plot(log(data_ss(:, yCol)), log(data_ss(:,ameanCol)),'.');
set(gca, 'FontSize', 14);
ylabel('log(Amean) (ms-2)', 'FontSize', 14);
xlabel('log(y) (m)', 'FontSize', 14);
subplot(1,2,2)
plot((data_ss(:, yCol)), (data_ss(:,ameanCol)),'.');
set(gca, 'FontSize', 14);
ylabel('Amean (ms-2)', 'FontSize', 14);
xlabel('y (m)', 'FontSize', 14);

figure;
subplot(1,2,1)
plot(log(data_ss(:, deltarCol)), log(data_ss(:,ameanCol)),'.');
set(gca, 'FontSize', 14);
ylabel('log(Amean) (ms-2)', 'FontSize', 14);
xlabel('log(deltar) (s-1)', 'FontSize', 14);
subplot(1,2,2)
plot((data_ss(:, deltarCol)), (data_ss(:,ameanCol)),'.');
set(gca, 'FontSize', 14);
ylabel('Amean (ms-2)', 'FontSize', 14);
xlabel('deltar (s-1)', 'FontSize', 14);

figure;
subplot(1,2,1)
plot(log(data_ss(:, rrefCol)), log(data_ss(:,ameanCol)),'.');
set(gca, 'FontSize', 14);
ylabel('log(Amean) (ms-2)', 'FontSize', 14);
xlabel('log(rref) (s-1)', 'FontSize', 14);
subplot(1,2,2)
plot((data_ss(:, rrefCol)), (data_ss(:,ameanCol)),'.');
set(gca, 'FontSize', 14);
ylabel('Amean (ms-2)', 'FontSize', 14);
xlabel('rref (s-1)', 'FontSize', 14);

% figure;
% subplot(1,2,1)
% plot(log(data_ss(:, r0Col)), log(data_ss(:,rdotCol)),'.');
% set(gca, 'FontSize', 14);
% ylabel('log(rdot) (s-2)', 'FontSize', 14);
% xlabel('log(r0) (s-1)', 'FontSize', 14);
% subplot(1,2,2)
% plot((data_ss(:, r0Col)), (data_ss(:,rdotCol)),'.');
% set(gca, 'FontSize', 14);
% ylabel('rdot (s-2)', 'FontSize', 14);
% xlabel('r0 (s-1)', 'FontSize', 14);

%% Plot Amean as a function of y, delta_r and "fixed" rref (FINAL - TO BE USED)
% One panel for the main figures (y on continuous scale)

% deltar as x1, y as x2 and rref as x3
close all;
rdotCol = ameanCol; % Shortcut (correct this later)
rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); wind = data_ss(:,windCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);

% values for which lines will be plotted
% y_percentile_5_50_95 = quantile(y, [0.05, 0.5, 0.95]);
% rref_percentile_5_50_95 = quantile(rref, [0.05, 0.5, 0.95]);
% 
% Untransformed domain
res = rdot; x1 = deltar; x2 = y; x3 = rref; x4 = wind;

x3_centers = quantile(x3,[0.15 0.5 0.85]);% quantile(x3,[0.25 0.5 0.75]); % values for which rref data is binned
binwidth = 1; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)

x2_forplot = [];
for ct=1:3 % for each wind condition
    for ct1=1:length(x3_centers)  % for each x3
        dummy = x2(data_ss(:,windCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2));
        x2_forplot = [x2_forplot; dummy]; 
    end
end
% x2range = [min(x2_forplot) max(x2_forplot)];
% x2range = [round(x2range(1)*10)/10 round(max(x2range(2))*10)/10];
x2range = [0.05 0.45]; % y range

% cmap=brewermap(round(diff(x2range)/0.01),'Set1'); 
% cmap = jet(round(diff(x2range)/0.01));  
cmap = copper(round(diff(x2range)/0.01));  

% cmap = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880];
% % % % % % cmap = [253,208,162; 253,174,107; 253,141,60; 241,105,19; 217,72,1; 166,54,3]./255;
edges = round(linspace(x2range(1),x2range(2),size(cmap,1)+1),2); % # edges = # color bins + 1

% For plotting continuous lines
% x2quantiles = quantile(x2_forplot,[0.05 0.5 0.95]) 
x2quantiles = quantile(data_ss(:,yCol),[0.05 0.5 0.95]); % x2 values for which lines are plotted
x2quantiles = [0.1 0.2 0.3 0.4]; %quantile(data_ss(:,yCol),[0.05 0.5 0.95]) % x2 values for which lines are plotted
% deltar_values_cont = 0.2:0.01:4; % deltar values for which lines are plotted

figure2d = figure; hold on;
colormap(cmap);

figure2d_log = figure; hold on;
colormap(cmap);

for dum=3%1:3 % for each wind condition
    ct = chosen_winds(dum);
    for ct1=2%1:length(x3_centers)  % for each x3
        figure(figure2d); hold on;
%         subplot(3, length(x3_centers), (ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,windCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2),:);
        y_ss2 = data_ss2(:,yCol);
%         deltar_ss2 = data_ss2(:,deltarCol);
%         [min(deltar_ss2) max(deltar_ss2)]
        [data_cmap, ~] = discretize(y_ss2, edges);
        data_cmap(y_ss2<=edges(1)) = 1; data_cmap(y_ss2>=edges(end)) = size(cmap,1);
        scatter(data_ss2(:,deltarCol), data_ss2(:,rdotCol), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
%         x2quantiles = quantile(deltar_ss2,[0.05 0.5 0.95]) % deltar values for which lines are plotted
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        dum = [];
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(x2quantiles(ct2)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x3_centers(ct1)) + coeffs(7)*log(deltar_values_cont)*log(x2quantiles(ct2)) + ...
                + coeffs(8)*log(x3_centers(ct1))*log(deltar_values_cont);
            if ct == 4
                % coeffs = [intercept logy wind2 wind3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif ct == 6
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            dum(ct2) = linecolor;
            plot(deltar_values_cont, exp(rdot_response),'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        display(dum)
        
        set(gca, 'FontSize', 16);
%         ylim([0, 35]);
%         yticks([0:5:35]);
        if ct1==1
            xlim([0 2.5]);
            xticks([0:0.5:2.5]);
        elseif ct1==2
            xlim([0 3.3]);
            xticks([0:0.5:3]);
        elseif ct1==3
            xlim([0 4]);
            xticks([0:0.5:4]);
        end
%         figure2d.Position = [-1313 46 895 820];

        
        figure(figure2d_log); hold on;
%         subplot(3, length(x3_centers), (ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,windCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2),:);
        y_ss2 = data_ss2(:,yCol);
        [data_cmap, ~] = discretize(y_ss2, edges);
        data_cmap(y_ss2<=edges(1)) = 1; data_cmap(y_ss2>=edges(end)) = size(cmap,1);
        scatter(log(data_ss2(:,deltarCol)), log(data_ss2(:,rdotCol)), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(x2quantiles(ct2)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x3_centers(ct1)) + coeffs(7)*log(deltar_values_cont)*log(x2quantiles(ct2)) + ...
                + coeffs(8)*log(x3_centers(ct1))*log(deltar_values_cont);
            if ct == 4
                % coeffs = [intercept logy wind2 wind3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif ct == 6
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            plot(log(deltar_values_cont), rdot_response,'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        
        set(gca, 'FontSize', 16);
%         ylim([1, 4]);
%         yticks([1:1:4]);
        xlim([-1.5 1.5]);
        xticks([-1.5:0.5:1.5]);
%         figure2d_log.Position = [-1313 46 895 820];

    end
end
% figure(figure2d);
% cmap_bar = colorbar('eastoutside');
% caxis(edges([1 end]));
% cmap_bar.Ticks = [edges(1) edges(end)];
% cmap_bar.Label.String = 'delta_r';
% cmap_bar.Label.FontSize = 16;


figure;
colormap(cmap);
cmap_bar = colorbar('eastoutside');
caxis(edges([1 end]));
cmap_bar.Ticks = [edges(1) x2quantiles edges(end)];
cmap_bar.Label.String = 'y';
cmap_bar.Label.FontSize = 16;
ylabel('rdot (s-2)', 'FontSize', 14);
xlabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);


%% Plot rdot as a function of delta_r, rref and "fixed" y (TO BE USED)
% deltar as x1, rref as x2 and y as x3
close all;
rdotCol = ameanCol;
rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); wind = data_ss(:,windCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);

% values for which lines will be plotted
y_percentile_5_50_95 = quantile(y, [0.05, 0.5, 0.95]);
rref_percentile_5_50_95 = quantile(rref, [0.05, 0.5, 0.95]);

% Untransformed domain
res = rdot; x1 = deltar; x2 = rref; x3 = y; x4 = wind;

x3_centers = quantile(x3,[0.25 0.5 0.75]); % values for which rref data is binned
binwidth = 0.05; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)
% x3_centers = [0.1 0.2 0.3];
% binwidth = 0.1; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)

x2_forplot = [];
for dum=1:3 % for each wind condition
    ct = chosen_winds(dum);
    for ct1=1:length(x3_centers)  % for each x3
        dummy = x2(data_ss(:,windCol) == ct & (data_ss(:,yCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,yCol) <= x3_centers(ct1)+binwidth/2));
        x2_forplot = [x2_forplot; dummy]; 
    end
end
% x2range = rref_percentile_5_50_95([1 end]);
% x2range = [round(x2range(1)*10)/10 round(max(x2range(2))*10)/10];
x2range = [1 5.5];

% cmap=brewermap(round(diff(x2range)/0.5),'Set1'); 
% cmap = jet(round(diff(x2range)/0.1));  
% cmap = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880];
cmap = copper(round(diff(x2range)/0.15));  

edges = round(linspace(x2range(1),x2range(2),size(cmap,1)+1),2); % # edges = # color bins + 1

% For plotting continuous lines
% x2quantiles = quantile(x2_forplot,[0.05 0.5 0.95]) 
x2quantiles = quantile(data_ss(:,rrefCol),[0.05 0.5 0.95]) % x2 values for which lines are plotted
% deltar_values_cont = 0.2:0.01:4; % deltar values for which lines are plotted

figure2d = figure; hold on;
colormap(cmap);

figure2d_log = figure; hold on;
colormap(cmap);
for dummy_ct=1%1:3 % for each wind condition
    ct=chosen_winds(dummy_ct);
    for ct1=2%1:length(x3_centers)  % for each x3
        figure(figure2d);
%         subplot(3, length(x3_centers), (dummy_ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,windCol) == ct & (data_ss(:,yCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,yCol) <= x3_centers(ct1)+binwidth/2),:);
        rref_ss2 = data_ss2(:,rrefCol);
%         deltar_ss2 = data_ss2(:,deltarCol);
%         [min(deltar_ss2) max(deltar_ss2)]
        [data_cmap, ~] = discretize(rref_ss2, edges);
        data_cmap(rref_ss2<=edges(1)) = 1; data_cmap(rref_ss2>=edges(end)) = size(cmap,1);
        scatter(data_ss2(:,deltarCol), data_ss2(:,rdotCol), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
%         x2quantiles = quantile(deltar_ss2,[0.05 0.5 0.95]) % deltar values for which lines are plotted
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        dum = [];
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(x3_centers(ct1)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x2quantiles(ct2)) + coeffs(7)*log(deltar_values_cont)*log(x3_centers(ct1)) + ...
                + coeffs(8)*log(x2quantiles(ct2))*log(deltar_values_cont);
            if dummy_ct == 2
                % coeffs = [intercept logy wind2 wind3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif dummy_ct == 3
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            dum(ct2) = linecolor;
            plot(deltar_values_cont, exp(rdot_response),'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        display(dum)
        
        set(gca, 'FontSize', 16);
%         ylim([0, 35]);
%         yticks([0:5:35]);
        if ct1==1
            xlim([0 5]);
            xticks([0:1:5]);
        elseif ct1==2
            xlim([0 6]);
            xticks([0:1:6]);
        elseif ct1==3
            xlim([0 6]);
            xticks([0:1:6]);
        end
%         figure2d.Position = [-1313 46 895 820];
        
        figure(figure2d_log);
%         subplot(3, length(x3_centers), (dummy_ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,windCol) == ct & (data_ss(:,yCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,yCol) <= x3_centers(ct1)+binwidth/2),:);
        rref_ss2 = data_ss2(:,rrefCol);
        [data_cmap, ~] = discretize(rref_ss2, edges);
        data_cmap(rref_ss2<=edges(1)) = 1; data_cmap(rref_ss2>=edges(end)) = size(cmap,1);
        scatter(log(data_ss2(:,deltarCol)), log(data_ss2(:,rdotCol)), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(x3_centers(ct1)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x2quantiles(ct2)) + coeffs(7)*log(deltar_values_cont)*log(x3_centers(ct1)) + ...
                + coeffs(8)*log(x2quantiles(ct2))*log(deltar_values_cont);
            if dummy_ct == 2
                % coeffs = [intercept logy wind2 wind3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif dummy_ct == 3
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            plot(log(deltar_values_cont), rdot_response,'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        
        set(gca, 'FontSize', 16);
%         ylim([0, 4]);
%         yticks([0:1:4]);
        xlim([-1.5 2]);
        xticks([-1.5:0.5:2]);
%         figure2d_log.Position = [-1313 46 895 820];
    end
end
% figure(figure2d);
% cmap_bar = colorbar('eastoutside');
% caxis(edges([1 end]));
% cmap_bar.Ticks = [edges(1) edges(end)];
% cmap_bar.Label.String = 'delta_r';
% cmap_bar.Label.FontSize = 16;


figure;
colormap(cmap);
cmap_bar = colorbar('eastoutside');
caxis(edges([1 end]));
cmap_bar.Ticks = [edges(1) x2quantiles edges(end)];
% cmap_bar.Ticks = [edges(1) edges(end)];
cmap_bar.Label.String = 'rref';
cmap_bar.Label.FontSize = 16;
ylabel('rdot (s-2)', 'FontSize', 14);
xlabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);


%% Plot Amean as a function of y, delta_r and rref (FINAL - TO BE USED)
% deltar as x1, y as x2 and rref as x3
close all;
amean = data_ss(:,ameanCol); y = data_ss(:,yCol); light = data_ss(:,lightCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);

% values for which lines will be plotted
% y_percentile_5_50_95 = quantile(y, [0.05, 0.5, 0.95]);
% rref_percentile_5_50_95 = quantile(rref, [0.05, 0.5, 0.95]);
% 
% Untransformed domain
res = amean; x1 = deltar; x2 = y; x3 = rref; x4 = light;

x3_centers = quantile(x3,[0.15 0.5 0.85]);% quantile(x3,[0.25 0.5 0.75]); % values for which rref data is binned
binwidth = 1; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)

x2_forplot = [];
for ct=1:3 % for each light condition
    for ct1=1:length(x3_centers)  % for each x3
        dummy = x2(data_ss(:,lightCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2));
        x2_forplot = [x2_forplot; dummy]; 
    end
end
% x2range = [min(x2_forplot) max(x2_forplot)];
% x2range = [round(x2range(1)*10)/10 round(max(x2range(2))*10)/10];
x2range = [0.05 0.35]; % y range

% cmap=brewermap(round(diff(x2range)/0.5),'Set1'); 
% cmap = jet(round(diff(x2range)/0.1));  
cmap = copper(round(diff(x2range)/0.01));  

% cmap = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880];
% cmap = [253,208,162; 253,174,107; 253,141,60; 241,105,19; 217,72,1; 166,54,3]./255;
edges = round(linspace(x2range(1),x2range(2),size(cmap,1)+1),2); % # edges = # color bins + 1

% For plotting continuous lines
% x2quantiles = quantile(x2_forplot,[0.05 0.5 0.95]) 
x2quantiles = quantile(data_ss(:,yCol),[0.05 0.5 0.95]); % x2 values for which lines are plotted
x2quantiles = [0.1 0.2 0.3]; %quantile(data_ss(:,yCol),[0.05 0.5 0.95]) % x2 values for which lines are plotted
% deltar_values_cont = 0.2:0.01:4; % deltar values for which lines are plotted

figure2d = figure; hold on;
colormap(cmap);

figure2d_log = figure; hold on;
colormap(cmap);
for ct=3%1:3 % for each light condition
    for ct1=2%1:length(x3_centers)  % for each x3
        figure(figure2d); hold on;
%         subplot(3, length(x3_centers), (ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,lightCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2),:);
        y_ss2 = data_ss2(:,yCol);
%         deltar_ss2 = data_ss2(:,deltarCol);
%         [min(deltar_ss2) max(deltar_ss2)]
        [data_cmap, ~] = discretize(y_ss2, edges);
        data_cmap(y_ss2<=edges(1)) = 1; data_cmap(y_ss2>=edges(end)) = size(cmap,1);
        scatter(data_ss2(:,deltarCol), data_ss2(:,ameanCol), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
%         x2quantiles = quantile(deltar_ss2,[0.05 0.5 0.95]) % deltar values for which lines are plotted
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        dum = [];
        for ct2=1:length(x2quantiles)
            amean_response = coeffs(1) + coeffs(2)*log(x2quantiles(ct2)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x3_centers(ct1)) + coeffs(7)*log(deltar_values_cont)*log(x2quantiles(ct2)) + ...
                + coeffs(8)*log(x3_centers(ct1))*log(deltar_values_cont);
            if ct == 2
                % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                amean_response = amean_response + coeffs(3);
            elseif ct == 3
                amean_response = amean_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            dum(ct2) = linecolor;
            plot(deltar_values_cont, exp(amean_response),'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        display(dum)
        
        set(gca, 'FontSize', 16);
%         ylim([0, 35]);
%         yticks([0:5:35]);
%         if ct1==1
%             xlim([0 2.5]);
%             xticks([0:0.5:2.5]);
%         elseif ct1==2
            xlim([0.5 3]);
            xticks([0.5:0.5:3]);
%         elseif ct1==3
%             xlim([0 4]);
%             xticks([0:0.5:4]);
%         end
%         figure2d.Position = [-1313 46 895 820];

        
        figure(figure2d_log); hold on;
%         subplot(3, length(x3_centers), (ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,lightCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2),:);
        y_ss2 = data_ss2(:,yCol);
        [data_cmap, ~] = discretize(y_ss2, edges);
        data_cmap(y_ss2<=edges(1)) = 1; data_cmap(y_ss2>=edges(end)) = size(cmap,1);
        scatter(log(data_ss2(:,deltarCol)), log(data_ss2(:,ameanCol)), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        for ct2=1:length(x2quantiles)
            amean_response = coeffs(1) + coeffs(2)*log(x2quantiles(ct2)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x3_centers(ct1)) + coeffs(7)*log(deltar_values_cont)*log(x2quantiles(ct2)) + ...
                + coeffs(8)*log(x3_centers(ct1))*log(deltar_values_cont);
            if ct == 2
                % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                amean_response = amean_response + coeffs(3);
            elseif ct == 3
                amean_response = amean_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            plot(log(deltar_values_cont), amean_response,'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        
        set(gca, 'FontSize', 16);
%         ylim([1, 4]);
%         yticks([1:1:4]);
%         xlim([-1.5 1.5]);
%         xticks([-1.5:0.5:1.5]);
%         figure2d_log.Position = [-1313 46 895 820];

    end
end
% figure(figure2d);
% cmap_bar = colorbar('eastoutside');
% caxis(edges([1 end]));
% cmap_bar.Ticks = [edges(1) edges(end)];
% cmap_bar.Label.String = 'delta_r';
% cmap_bar.Label.FontSize = 16;


figure;
colormap(cmap);
cmap_bar = colorbar('eastoutside');
caxis(edges([1 end]));
cmap_bar.Ticks = [edges(1) x2quantiles edges(end)];
cmap_bar.Label.String = 'y';
cmap_bar.Label.FontSize = 16;
ylabel('Amean (ms-2)', 'FontSize', 14);
xlabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);

%% Plot Amean as a function of rref (FINAL TO BE USED)
% rref as x1
close all;
amean = data_ss(:,ameanCol); y = data_ss(:,yCol); light = data_ss(:,lightCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);

% values for which lines will be plotted
y_chosen = 0.21; %quantile(y, 0.5); % mean(y)
y_bin = 0.06; %mean(y) std(y)
deltar_chosen = 1.68; quantile(deltar, 0.5); %mean(delta_r) std(delta_r)
deltar_bin = 0.8;

% Untransformed domain
res = amean; x1 = rref;

% x3_centers = quantile(x3,[0.25 0.5 0.75]); % values for which rref data is binned
% binwidth = 0.05; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)

x1_forplot = [];
for ct=1:3 % for each light condition
    dummy = x1(data_ss(:,lightCol) == ct & (data_ss(:,yCol) >= y_chosen-y_bin/2 & data_ss(:,yCol) <= y_chosen+y_bin/2) & ...
                                           (data_ss(:,deltarCol) >= deltar_chosen-deltar_bin/2 & data_ss(:,deltarCol) <= deltar_chosen+deltar_bin/2));
    x1_forplot = [x1_forplot; dummy];
end
x1range = [min(x1_forplot) max(x1_forplot)];
x1range = [floor(x1range(1)*10)/10 round(max(x1range(2))*10)/10];
rref_values_cont = x1range(1):0.1:x1range(2);
% x2range = [1 5.5];

% cmap=brewermap(round(diff(x2range)/0.5),'Set1'); 
% cmap = jet(round(diff(x2range)/0.1));  
cmap = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880];
edges = round(linspace(x2range(1),x2range(2),size(cmap,1)+1),2); % # edges = # color bins + 1

% For plotting continuous lines
% x2quantiles = quantile(x2_forplot,[0.05 0.5 0.95]) 
% x2quantiles = quantile(data_ss(:,rrefCol),[0.05 0.5 0.95]) % x2 values for which lines are plotted
 % deltar values for which lines are plotted

figure2d = figure; hold on;
colormap(cmap);

figure2d_log = figure; hold on;
colormap(cmap);
for ct=3%1:3 % for each light condition
    
        figure(figure2d); hold on;
%         subplot(3, 1, ct); hold on;
        data_ss2 = data_ss(data_ss(:,lightCol) == ct & (data_ss(:,yCol) >= y_chosen-y_bin/2 & data_ss(:,yCol) <= y_chosen+y_bin/2) & ...
                                           (data_ss(:,deltarCol) >= deltar_chosen-deltar_bin/2 & data_ss(:,deltarCol) <= deltar_chosen+deltar_bin/2),:);
%         rref_ss2 = data_ss2(:,rrefCol);
%         [data_cmap, ~] = discretize(rref_ss2, edges);
%         data_cmap(rref_ss2<=edges(1)) = 1; data_cmap(rref_ss2>=edges(end)) = size(cmap,1);
        scatter(data_ss2(:,rrefCol), data_ss2(:,ameanCol), 10,'filled','o');
        
        % Plot continous lines
        
        amean_response = coeffs(1) + coeffs(2)*log(y_chosen) + coeffs(5)*log(deltar_chosen) + ...
            coeffs(6)*log(rref_values_cont) + coeffs(7)*log(y_chosen)*log(deltar_chosen) + ...
            + coeffs(8)*log(deltar_chosen)*log(rref_values_cont);
        if ct == 2
            % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
            amean_response = amean_response + coeffs(3);
        elseif ct == 3
            amean_response = amean_response + coeffs(4);
        end
%         linecolor = discretize(x2quantiles(ct2), edges);
        plot(rref_values_cont, exp(amean_response),'Color',[0 0 0],'LineWidth', 2);
        
        set(gca, 'FontSize', 16);
        ylim([0, 4]);
        yticks([0:1:4]);
%         if ct1==1
%             xlim([1.5 5]);
%             xticks([1:1:5]);
%         elseif ct1==2
%             xlim([0 6]);
%             xticks([0:1:6]);
%         elseif ct1==3
%             xlim([0 6]);
%             xticks([0:1:6]);
%         end
%         figure2d.Position = [326    42   335   954];
        
        figure(figure2d_log); hold on;
%         subplot(3, 1, ct); hold on;
        data_ss2 = data_ss(data_ss(:,lightCol) == ct & (data_ss(:,yCol) >= y_chosen-y_bin/2 & data_ss(:,yCol) <= y_chosen+y_bin/2) & ...
                                           (data_ss(:,deltarCol) >= deltar_chosen-deltar_bin/2 & data_ss(:,deltarCol) <= deltar_chosen+deltar_bin/2),:);
        scatter(log(data_ss2(:,rrefCol)), log(data_ss2(:,ameanCol)), 10,'filled','o');
        
        % Plot continous lines
        plot(log(rref_values_cont), amean_response,'Color',[0 0 0],'LineWidth', 2);
        
        
        set(gca, 'FontSize', 16);
%         ylim([1.5 3]);
%         yticks([1.5:0.5:3]);
%         xlim([0.5 1.75]);
%         xticks([0.5:0.5:1.5]);
%         figure2d_log.Position = [667    34   335   961];
end
% figure(figure2d);
% cmap_bar = colorbar('eastoutside');
% caxis(edges([1 end]));
% cmap_bar.Ticks = [edges(1) edges(end)];
% cmap_bar.Label.String = 'delta_r';
% cmap_bar.Label.FontSize = 16;


figure;
% colormap(cmap);
% cmap_bar = colorbar('eastoutside');
% caxis(edges([1 end]));
% cmap_bar.Ticks = [edges(1) edges(end)];
% cmap_bar.Label.String = 'rref';
% cmap_bar.Label.FontSize = 16;
ylabel('rdot (s-2)', 'FontSize', 14);
xlabel('rref (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);

%% Plot Amean as a function of y, delta_r and rref (FINAL - TO BE USED) - supplement file
% deltar as x1, y as x2 and rref as x3
close all;
amean = data_ss(:,ameanCol); y = data_ss(:,yCol); light = data_ss(:,lightCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);

% values for which lines will be plotted
% y_percentile_5_50_95 = quantile(y, [0.05, 0.5, 0.95]);
% rref_percentile_5_50_95 = quantile(rref, [0.05, 0.5, 0.95]);
% 
% Untransformed domain
res = amean; x1 = deltar; x2 = y; x3 = rref; x4 = light;

x3_centers = quantile(x3,[0.15 0.5 0.85]);% quantile(x3,[0.25 0.5 0.75]); % values for which rref data is binned
binwidth = 1; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)

x2_forplot = [];
for ct=1:3 % for each light condition
    for ct1=1:length(x3_centers)  % for each x3
        dummy = x2(data_ss(:,lightCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2));
        x2_forplot = [x2_forplot; dummy]; 
    end
end
% x2range = [min(x2_forplot) max(x2_forplot)];
% x2range = [round(x2range(1)*10)/10 round(max(x2range(2))*10)/10];
x2range = [0.05 0.35]; % y range

% cmap=brewermap(round(diff(x2range)/0.5),'Set1'); 
% cmap = jet(round(diff(x2range)/0.1));  
cmap = copper(round(diff(x2range)/0.01));  

% cmap = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880];
% cmap = [253,208,162; 253,174,107; 253,141,60; 241,105,19; 217,72,1; 166,54,3]./255;
edges = round(linspace(x2range(1),x2range(2),size(cmap,1)+1),2); % # edges = # color bins + 1

% For plotting continuous lines
% x2quantiles = quantile(x2_forplot,[0.05 0.5 0.95]) 
x2quantiles = quantile(data_ss(:,yCol),[0.05 0.5 0.95]); % x2 values for which lines are plotted
x2quantiles = [0.1 0.2 0.3]; %quantile(data_ss(:,yCol),[0.05 0.5 0.95]) % x2 values for which lines are plotted
% deltar_values_cont = 0.2:0.01:4; % deltar values for which lines are plotted

figure2d = figure; hold on;
colormap(cmap);

figure2d_log = figure; hold on;
colormap(cmap);
for ct=1:3 % for each light condition
    for ct1=1:length(x3_centers)  % for each x3
        figure(figure2d); hold on;
        subplot(3, length(x3_centers), (ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,lightCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2),:);
        y_ss2 = data_ss2(:,yCol);
%         deltar_ss2 = data_ss2(:,deltarCol);
%         [min(deltar_ss2) max(deltar_ss2)]
        [data_cmap, ~] = discretize(y_ss2, edges);
        data_cmap(y_ss2<=edges(1)) = 1; data_cmap(y_ss2>=edges(end)) = size(cmap,1);
        scatter(data_ss2(:,deltarCol), data_ss2(:,ameanCol), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
%         x2quantiles = quantile(deltar_ss2,[0.05 0.5 0.95]) % deltar values for which lines are plotted
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        dum = [];
        for ct2=1:length(x2quantiles)
            amean_response = coeffs(1) + coeffs(2)*log(x2quantiles(ct2)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x3_centers(ct1)) + coeffs(7)*log(deltar_values_cont)*log(x2quantiles(ct2)) + ...
                + coeffs(8)*log(x3_centers(ct1))*log(deltar_values_cont);
            if ct == 2
                % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                amean_response = amean_response + coeffs(3);
            elseif ct == 3
                amean_response = amean_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            dum(ct2) = linecolor;
            plot(deltar_values_cont, exp(amean_response),'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        display(dum)
        
        set(gca, 'FontSize', 16);
%         ylim([0, 35]);
%         yticks([0:5:35]);
%         if ct1==1
%             xlim([0 2.5]);
%             xticks([0:0.5:2.5]);
%         elseif ct1==2
            xlim([0.5 3]);
            xticks([0.5:0.5:3]);
%         elseif ct1==3
%             xlim([0 4]);
%             xticks([0:0.5:4]);
%         end
%         figure2d.Position = [-1313 46 895 820];

        
        figure(figure2d_log); hold on;
        subplot(3, length(x3_centers), (ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,lightCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2),:);
        y_ss2 = data_ss2(:,yCol);
        [data_cmap, ~] = discretize(y_ss2, edges);
        data_cmap(y_ss2<=edges(1)) = 1; data_cmap(y_ss2>=edges(end)) = size(cmap,1);
        scatter(log(data_ss2(:,deltarCol)), log(data_ss2(:,ameanCol)), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        for ct2=1:length(x2quantiles)
            amean_response = coeffs(1) + coeffs(2)*log(x2quantiles(ct2)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x3_centers(ct1)) + coeffs(7)*log(deltar_values_cont)*log(x2quantiles(ct2)) + ...
                + coeffs(8)*log(x3_centers(ct1))*log(deltar_values_cont);
            if ct == 2
                % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                amean_response = amean_response + coeffs(3);
            elseif ct == 3
                amean_response = amean_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            plot(log(deltar_values_cont), amean_response,'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        
        set(gca, 'FontSize', 16);
%         ylim([1, 4]);
%         yticks([1:1:4]);
%         xlim([-1.5 1.5]);
%         xticks([-1.5:0.5:1.5]);
%         figure2d_log.Position = [-1313 46 895 820];

    end
end
% figure(figure2d);
% cmap_bar = colorbar('eastoutside');
% caxis(edges([1 end]));
% cmap_bar.Ticks = [edges(1) edges(end)];
% cmap_bar.Label.String = 'delta_r';
% cmap_bar.Label.FontSize = 16;


figure;
colormap(cmap);
cmap_bar = colorbar('eastoutside');
caxis(edges([1 end]));
cmap_bar.Ticks = [edges(1) x2quantiles edges(end)];
cmap_bar.Label.String = 'y';
cmap_bar.Label.FontSize = 16;
ylabel('Amean (ms-2)', 'FontSize', 14);
xlabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);

%% Plot Amean as a function of rref (FINAL TO BE USED) - supplement file
% rref as x1
close all;
amean = data_ss(:,ameanCol); y = data_ss(:,yCol); light = data_ss(:,lightCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);

% values for which lines will be plotted
y_chosen = 0.21; %quantile(y, 0.5); % mean(y)
y_bin = 0.06; %mean(y) std(y)
deltar_chosen = 1.68; quantile(deltar, 0.5); %mean(delta_r) std(delta_r)
deltar_bin = 0.8;

% Untransformed domain
res = amean; x1 = rref;

% x3_centers = quantile(x3,[0.25 0.5 0.75]); % values for which rref data is binned
% binwidth = 0.05; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)

x1_forplot = [];
for ct=1:3 % for each light condition
    dummy = x1(data_ss(:,lightCol) == ct & (data_ss(:,yCol) >= y_chosen-y_bin/2 & data_ss(:,yCol) <= y_chosen+y_bin/2) & ...
                                           (data_ss(:,deltarCol) >= deltar_chosen-deltar_bin/2 & data_ss(:,deltarCol) <= deltar_chosen+deltar_bin/2));
    x1_forplot = [x1_forplot; dummy];
end
x1range = [min(x1_forplot) max(x1_forplot)];
x1range = [floor(x1range(1)*10)/10 round(max(x1range(2))*10)/10];
rref_values_cont = x1range(1):0.1:x1range(2);
% x2range = [1 5.5];

% cmap=brewermap(round(diff(x2range)/0.5),'Set1'); 
% cmap = jet(round(diff(x2range)/0.1));  
cmap = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880];
edges = round(linspace(x2range(1),x2range(2),size(cmap,1)+1),2); % # edges = # color bins + 1

% For plotting continuous lines
% x2quantiles = quantile(x2_forplot,[0.05 0.5 0.95]) 
% x2quantiles = quantile(data_ss(:,rrefCol),[0.05 0.5 0.95]) % x2 values for which lines are plotted
 % deltar values for which lines are plotted

figure2d = figure; hold on;
colormap(cmap);

figure2d_log = figure; hold on;
colormap(cmap);
for ct=1:3 % for each light condition
    
        figure(figure2d); hold on;
        subplot(3, 1, ct); hold on;
        data_ss2 = data_ss(data_ss(:,lightCol) == ct & (data_ss(:,yCol) >= y_chosen-y_bin/2 & data_ss(:,yCol) <= y_chosen+y_bin/2) & ...
                                           (data_ss(:,deltarCol) >= deltar_chosen-deltar_bin/2 & data_ss(:,deltarCol) <= deltar_chosen+deltar_bin/2),:);
%         rref_ss2 = data_ss2(:,rrefCol);
%         [data_cmap, ~] = discretize(rref_ss2, edges);
%         data_cmap(rref_ss2<=edges(1)) = 1; data_cmap(rref_ss2>=edges(end)) = size(cmap,1);
        scatter(data_ss2(:,rrefCol), data_ss2(:,ameanCol), 10,'filled','o');
        
        % Plot continous lines
        
        amean_response = coeffs(1) + coeffs(2)*log(y_chosen) + coeffs(5)*log(deltar_chosen) + ...
            coeffs(6)*log(rref_values_cont) + coeffs(7)*log(y_chosen)*log(deltar_chosen) + ...
            + coeffs(8)*log(deltar_chosen)*log(rref_values_cont);
        if ct == 2
            % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
            amean_response = amean_response + coeffs(3);
        elseif ct == 3
            amean_response = amean_response + coeffs(4);
        end
%         linecolor = discretize(x2quantiles(ct2), edges);
        plot(rref_values_cont, exp(amean_response),'Color',[0 0 0],'LineWidth', 2);
        
        set(gca, 'FontSize', 16);
        ylim([0, 4]);
        yticks([0:1:4]);
%         if ct1==1
%             xlim([1.5 5]);
%             xticks([1:1:5]);
%         elseif ct1==2
%             xlim([0 6]);
%             xticks([0:1:6]);
%         elseif ct1==3
%             xlim([0 6]);
%             xticks([0:1:6]);
%         end
%         figure2d.Position = [326    42   335   954];
        
        figure(figure2d_log); hold on;
        subplot(3, 1, ct); hold on;
        data_ss2 = data_ss(data_ss(:,lightCol) == ct & (data_ss(:,yCol) >= y_chosen-y_bin/2 & data_ss(:,yCol) <= y_chosen+y_bin/2) & ...
                                           (data_ss(:,deltarCol) >= deltar_chosen-deltar_bin/2 & data_ss(:,deltarCol) <= deltar_chosen+deltar_bin/2),:);
        scatter(log(data_ss2(:,rrefCol)), log(data_ss2(:,ameanCol)), 10,'filled','o');
        
        % Plot continous lines
        plot(log(rref_values_cont), amean_response,'Color',[0 0 0],'LineWidth', 2);
        
        
        set(gca, 'FontSize', 16);
%         ylim([1.5 3]);
%         yticks([1.5:0.5:3]);
%         xlim([0.5 1.75]);
%         xticks([0.5:0.5:1.5]);
%         figure2d_log.Position = [667    34   335   961];
end
% figure(figure2d);
% cmap_bar = colorbar('eastoutside');
% caxis(edges([1 end]));
% cmap_bar.Ticks = [edges(1) edges(end)];
% cmap_bar.Label.String = 'delta_r';
% cmap_bar.Label.FontSize = 16;


figure;
% colormap(cmap);
% cmap_bar = colorbar('eastoutside');
% caxis(edges([1 end]));
% cmap_bar.Ticks = [edges(1) edges(end)];
% cmap_bar.Label.String = 'rref';
% cmap_bar.Label.FontSize = 16;
ylabel('rdot (s-2)', 'FontSize', 14);
xlabel('rref (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);


%% Plot rdot as a function of y, delta_r and rref (TO BE USED)
% deltar as x1, rref as x2 and y as x3
close all;
rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); light = data_ss(:,lightCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);

% values for which lines will be plotted
y_percentile_5_50_95 = quantile(y, [0.05, 0.5, 0.95]);
rref_percentile_5_50_95 = quantile(rref, [0.05, 0.5, 0.95]);

% Untransformed domain
res = rdot; x1 = deltar; x2 = rref; x3 = y; x4 = light;

x3_centers = quantile(x3,[0.25 0.5 0.75]); % values for which rref data is binned
binwidth = 0.05; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)
% x3_centers = [0.1 0.2 0.3];
% binwidth = 0.1; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)

x2_forplot = [];
for ct=1:3 % for each light condition
    for ct1=1:length(x3_centers)  % for each x3
        dummy = x2(data_ss(:,lightCol) == ct & (data_ss(:,yCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,yCol) <= x3_centers(ct1)+binwidth/2));
        x2_forplot = [x2_forplot; dummy]; 
    end
end
% x2range = rref_percentile_5_50_95([1 end]);
% x2range = [round(x2range(1)*10)/10 round(max(x2range(2))*10)/10];
x2range = [1 5.5];

% cmap=brewermap(round(diff(x2range)/0.5),'Set1'); 
% cmap = jet(round(diff(x2range)/0.1));  
cmap = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880];
edges = round(linspace(x2range(1),x2range(2),size(cmap,1)+1),2); % # edges = # color bins + 1

% For plotting continuous lines
% x2quantiles = quantile(x2_forplot,[0.05 0.5 0.95]) 
x2quantiles = quantile(data_ss(:,rrefCol),[0.05 0.5 0.95]) % x2 values for which lines are plotted
% deltar_values_cont = 0.2:0.01:4; % deltar values for which lines are plotted

figure2d = figure; hold on;
colormap(cmap);

figure2d_log = figure; hold on;
colormap(cmap);
for ct=1:3 % for each light condition
    for ct1=1:length(x3_centers)  % for each x3
        figure(figure2d);
        subplot(3, length(x3_centers), (ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,lightCol) == ct & (data_ss(:,yCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,yCol) <= x3_centers(ct1)+binwidth/2),:);
        rref_ss2 = data_ss2(:,rrefCol);
%         deltar_ss2 = data_ss2(:,deltarCol);
%         [min(deltar_ss2) max(deltar_ss2)]
        [data_cmap, ~] = discretize(rref_ss2, edges);
        data_cmap(rref_ss2<=edges(1)) = 1; data_cmap(rref_ss2>=edges(end)) = size(cmap,1);
        scatter(data_ss2(:,deltarCol), data_ss2(:,rdotCol), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
%         x2quantiles = quantile(deltar_ss2,[0.05 0.5 0.95]) % deltar values for which lines are plotted
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        dum = [];
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(x3_centers(ct1)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x2quantiles(ct2)) + coeffs(7)*log(deltar_values_cont)*log(x3_centers(ct1)) + ...
                + coeffs(8)*log(x2quantiles(ct2))*log(deltar_values_cont);
            if ct == 2
                % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif ct == 3
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            dum(ct2) = linecolor;
            plot(deltar_values_cont, exp(rdot_response),'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        display(dum)
        
        set(gca, 'FontSize', 16);
        ylim([0, 35]);
        yticks([0:5:35]);
        if ct1==1
            xlim([0 5]);
            xticks([0:1:5]);
        elseif ct1==2
            xlim([0 6]);
            xticks([0:1:6]);
        elseif ct1==3
            xlim([0 6]);
            xticks([0:1:6]);
        end
        figure2d.Position = [-1313 46 895 820];
        
        figure(figure2d_log);
        subplot(3, length(x3_centers), (ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,lightCol) == ct & (data_ss(:,yCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,yCol) <= x3_centers(ct1)+binwidth/2),:);
        rref_ss2 = data_ss2(:,rrefCol);
        [data_cmap, ~] = discretize(rref_ss2, edges);
        data_cmap(rref_ss2<=edges(1)) = 1; data_cmap(rref_ss2>=edges(end)) = size(cmap,1);
        scatter(log(data_ss2(:,deltarCol)), log(data_ss2(:,rdotCol)), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(x3_centers(ct1)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x2quantiles(ct2)) + coeffs(7)*log(deltar_values_cont)*log(x3_centers(ct1)) + ...
                + coeffs(8)*log(x2quantiles(ct2))*log(deltar_values_cont);
            if ct == 2
                % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif ct == 3
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            plot(log(deltar_values_cont), rdot_response,'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        
        set(gca, 'FontSize', 16);
        ylim([0, 4]);
        yticks([0:1:4]);
        xlim([-1.5 2]);
        xticks([-1.5:0.5:2]);
        figure2d_log.Position = [-1313 46 895 820];
    end
end
% figure(figure2d);
% cmap_bar = colorbar('eastoutside');
% caxis(edges([1 end]));
% cmap_bar.Ticks = [edges(1) edges(end)];
% cmap_bar.Label.String = 'delta_r';
% cmap_bar.Label.FontSize = 16;


figure;
colormap(cmap);
cmap_bar = colorbar('eastoutside');
caxis(edges([1 end]));
cmap_bar.Ticks = [edges(1) edges(end)];
cmap_bar.Label.String = 'rref';
cmap_bar.Label.FontSize = 16;
ylabel('rdot (s-2)', 'FontSize', 14);
xlabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);