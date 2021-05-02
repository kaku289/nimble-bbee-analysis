%% Plot results from Amean model
%% Plot statistical model from R
close all;
% R model
% Factor = 1.5, with ymean as y
% coeffs = [1.41786160 -0.23188443  0.03265318  0.10269169 -0.12191186  0.30937950 -0.36988825 -0.11229241]; % [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]

% Factor = 1.5, with ystart as y
% coeffs = [1.33920177 -0.29558317  0.03165654  0.09682028 -0.15971722  0.30151310 -0.35868115 0]; % [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]

% Factor = 1, with ystart as y
coeffs = [1.58702290  0.42800399  0.05374526  0.16907431  0.32467186 -1.45318255 -0.58756146 0];
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
patternCol = 3;
lightCol = 4;
timeCol = 5;
dayCol = 6;
yCol = 14;
rrefCol = 8;
rdotCol = 9;
factorCol = 10;
isriseCol = 11;
deltarCol = 12;
deltaVCol = 15;
deltatCol = 16;
ameanCol = 17;
% data_write(:,end+1) = data_write(:,rrefCol)-data_write(:,deltarCol);
% r0Col = size(data_write,2);
data_write(:,end+1) = data_write(:,deltaVCol)./data_write(:,deltatCol);
dvdtCol = size(data_write,2);
% rdotCol = ameanCol; %dvdtCol;
ameanCol = 18;

% get subset of data

chosen_fac = 1;
% data_ss = data_write(abs(data_write(:,factorCol) - chosen_fac) < 1e-6 & data_write(:,dayCol) <= 20190710, :);
% data_ss = data_write(data_write(:,isriseCol) == 1 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6 & data_write(:,dayCol) <= 20190710, :);
data_ss = data_write(data_write(:,isriseCol) == 1 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6 & data_write(:,dayCol) <= 20190710 & data_write(:,ameanCol) > 0, :);

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

%% Plot histogram of Amean
chosen_fac = 1;
dayCol = 6;
rdotCol =9;
factorCol = 10;
isriseCol = 11;
ameanCol = 18;
data_ss = data_write(data_write(:,isriseCol) == 1 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6 & data_write(:,dayCol) <= 20190710, :);

figure;
histogram(data_ss(:,ameanCol));
% histfit(data_ss(:,ameanCol),[],'Normal')
% xlim([0 40]);
% histogram(-vertcat(data.rmean), [0:0.5:8]);
% histogram(-vertcat(data.rmean), [0:0.5:9.5]);

xlabel('Estimated rate-of-relative-rate-of-expansion, rdote (s-2)', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);
pd = fitdist(data_ss(:,rdotCol),'Gamma');
median(data_ss(:,rdotCol))

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
cmap = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880];
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
y_chosen = 0.22; %quantile(y, 0.5); % mean(y)
y_bin = 0.06; %mean(y) std(y)
deltar_chosen = 1.67; quantile(deltar, 0.5); %mean(delta_r) std(delta_r)
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

%% Mean accleration plot (FINAL - TO BE USED) - with rref and y on two axes and for a fixed deltar
% close all;
ameanCol = 17;
ameanCol = 18;

rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); light = data_ss(:,lightCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);
chosen_deltar = 1;%[1.6 2.4]%quantile(deltar,[0.25 0.5 0.75]);
rref_grid = 1:0.25:5;%quantile(rref,[0.05 0.95])
y_grid = 0.1:0.025:0.35;%quantile(y,[0.05 0.95])
rdot_response1 = nan(length(y_grid),length(rref_grid));
Amean_response2 = nan(length(y_grid),length(rref_grid)); % from rdot model
Amean_predicted = nan(length(y_grid),length(rref_grid)); % from Amean model
% coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
% coeffs_Amean = [1.71368183  0.46295979  0.04482469  0.14650472  0.44938972 -1.53883634 -0.55715188 0];


for ct_y=1:length(y_grid)
    for ct_rref=1:length(rref_grid)
        y0 = y_grid(ct_y); t0 = 0; thisRref = rref_grid(ct_rref);
        thisDeltar = chosen_deltar;
        r0 = thisRref - chosen_deltar;
        if r0 > 0
            V0 = r0*y0;
            tspan = [0:0.01:1.5];
            ystart = [y0 V0]; 
            % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
            log_rdot = coeffs(1) + coeffs(2)*log(y0) + coeffs(4) + coeffs(5)*log(thisDeltar) + ...
                coeffs(6)*log(thisRref) + coeffs(7)*log(y0)*log(thisDeltar) + ...
                + coeffs(8)*log(thisRref)*log(thisDeltar);
            rdot = exp(log_rdot);
            [t,y] = ode45(@(t,y) filteredState_BlindLandingtrack.const_drdt_movement(t,y,rdot), tspan, ystart);
            r = y(:,2)./y(:,1);
            tv_interp = interp1(r,[t y],thisRref);
            if tv_interp(2) > 0 && r(end) >= thisRref
                deltat = tv_interp(1); deltaV = tv_interp(3)-V0;
                Amean_response2(ct_y, ct_rref) = deltaV/deltat;
                
                rdot_response1(ct_y, ct_rref) = rdot;
                
%                 log_Amean_pred = coeffs_Amean(1) + coeffs_Amean(2)*log(y0) + coeffs_Amean(4) + coeffs_Amean(5)*log(thisDeltar) + ...
%                 coeffs_Amean(6)*log(chosen_rref) + coeffs_Amean(7)*log(y0)*log(thisDeltar) + ...
%                 + coeffs_Amean(8)*log(chosen_rref)*log(thisDeltar);
%                 Amean_predicted(ct_y, ct_deltar) = exp(log_Amean_pred);
            end
            
            
            
            
        else
            continue;
        end
    end
end

% Collect data points for plotting
binwidth = 0.2; 
data_ss2 = data_ss(data_ss(:,lightCol) == 3 & (data_ss(:,deltarCol) >= chosen_deltar-binwidth/2 & data_ss(:,deltarCol) <= chosen_deltar+binwidth/2),:);
% data_ss2 = data_ss((data_ss(:,deltarCol) >= chosen_deltar-binwidth/2 & data_ss(:,deltarCol) <= chosen_deltar+binwidth/2),:);
rdot = data_ss2(:,rdotCol); y = data_ss2(:,yCol); light = data_ss2(:,lightCol); deltar = data_ss2(:,deltarCol); rref = data_ss2(:,rrefCol); amean = data_ss2(:,ameanCol);

% managing colorbar
% rdot_range = [8 25];
rdot_range = [5 20];
nColors = round(diff(rdot_range)/1);
rdot_cmap = jet(2*nColors);
rdot_cmap = rdot_cmap(round(nColors/2):round(nColors/2)+nColors-1,:);
rdot_edges = round(linspace(rdot_range(1),rdot_range(2),size(rdot_cmap,1)+1),2); % # edges = # color bins + 1

% amean_range = [0 10];
amean_range = [0 2.5];
nColors = round(diff(amean_range)/0.5);
amean_cmap = jet(2*nColors);
amean_cmap = amean_cmap(round(nColors/2):round(nColors/2)+nColors-1,:);
% amean_cmap = jet(round(diff(amean_range)/0.5));  
amean_edges = round(linspace(amean_range(1),amean_range(2),size(amean_cmap,1)+1),2); % # edges = # color bins + 1

[Y, RREF] = meshgrid(rref_grid, y_grid);
rdot_fig = figure; hold on;
% subplot(1,2,1);
colormap(rdot_cmap);
[~,c1] = contourf(RREF, Y, rdot_response1, rdot_edges, 'LineStyle', 'none');

[data_cmap, ~] = discretize(rdot, rdot_edges);
data_cmap(rdot<=rdot_edges(1)) = 1; data_cmap(rdot>=rdot_edges(end)) = size(rdot_cmap,1);
scatter(y, rref, 25, rdot_cmap(data_cmap',:),'filled','o','MarkerEdgeColor',[0 .5 .5]);

title('rdot (s-2)', 'FontSize', 14);
xlabel('y (m)', 'FontSize', 14);
ylabel('rref (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);
xlim([0.1 0.35]);
xticks([0.1:0.05:0.35]);
ylim([2 5]);
yticks([2:0.5:5]);

cmap_bar = colorbar('eastoutside');
caxis(rdot_edges([1 end]));
cmap_bar.Ticks = [rdot_edges(1) rdot_edges(end)];
% cmap_bar.Label.String = 'rdot';
cmap_bar.Label.FontSize = 16;

amean_fig = figure; hold on;
% subplot(1,2,2) % rdot=f(y, deltar)
colormap(amean_cmap);
contourf(RREF, Y, Amean_response2, amean_edges, 'LineStyle', 'none');

[data_cmap, ~] = discretize(amean, amean_edges);
data_cmap(amean<=amean_edges(1)) = 1; data_cmap(amean>=amean_edges(end)) = size(amean_cmap,1);
scatter(y, rref, 25, amean_cmap(data_cmap',:),'filled','o','MarkerEdgeColor',[0 .5 .5]);

title('Amean (ms-2)', 'FontSize', 14);
xlabel('y (m)', 'FontSize', 14);
ylabel('rref (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);
xlim([0.1 0.35]);
xticks([0.1:0.05:0.35]);
ylim([2 5]);
yticks([2:0.5:5]);

cmap_bar = colorbar('eastoutside');
caxis(amean_edges([1 end]));
cmap_bar.Ticks = [amean_edges(1) amean_edges(end)];
% cmap_bar.Label.String = 'rdot';
cmap_bar.Label.FontSize = 16;

%% Mean accleration plot (3D data points - TO BE USED) - with rref, y and delta_r on three axes
close all;
ameanCol = 17;
ameanCol = 18;

% rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); light = data_ss(:,lightCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);

% Collect data points for plotting
binwidth = 0.2; 
data_ss2 = data_ss(data_ss(:,lightCol) == 3 & data_ss(:,isriseCol) == 1,:);
rdot = data_ss2(:,rdotCol); y = data_ss2(:,yCol); light = data_ss2(:,lightCol); deltar = data_ss2(:,deltarCol); rref = data_ss2(:,rrefCol); amean = data_ss2(:,ameanCol);

% quantile(amean,[0.05, 0.5, 0.95])
% figure;
% histogram(amean);

amean_range = [0 2.5];
nColors = round(diff(amean_range)/0.5);
amean_cmap = jet(2*nColors);
amean_cmap = amean_cmap(round(nColors/2):round(nColors/2)+nColors-1,:);
% amean_cmap = jet(round(diff(amean_range)/0.5));  
amean_edges = round(linspace(amean_range(1),amean_range(2),size(amean_cmap,1)+1),2); % # edges = # color bins + 1

amean_fig = figure; hold on;
% subplot(1,2,2) % rdot=f(y, deltar)
colormap(amean_cmap);

[data_cmap, ~] = discretize(amean, amean_edges);
data_cmap(amean<=amean_edges(1)) = 1; data_cmap(amean>=amean_edges(end)) = size(amean_cmap,1);
scatter3(y, rref, deltar, 25, amean_cmap(data_cmap',:),'filled','o','MarkerEdgeColor',[0 .5 .5]);

% title('Amean (ms-2)', 'FontSize', 14);
xlabel('y (m)', 'FontSize', 14);
ylabel('rref (s-1)', 'FontSize', 14);
zlabel('delta r (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);
% xlim([0.1 0.35]);
% xticks([0.1:0.05:0.35]);
% ylim([2 5]);
% yticks([2:0.5:5]);

cmap_bar = colorbar('eastoutside');
caxis(amean_edges([1 end]));
cmap_bar.Ticks = [amean_edges(1) amean_edges(end)];
cmap_bar.Label.String = 'Amean (ms-2)';
cmap_bar.Label.FontSize = 16;

view([90 0])
% figure;
% scatter(log(y), log(rref), 25, amean_cmap(data_cmap',:),'filled','o','MarkerEdgeColor',[0 .5 .5]);
% ylim([-3 4]);

%% Mean accleration plot (TO BE USED) - with deltar and y on two axes and for a fixed rref
close all;
% ameanCol = 17;
ameanCol = 18;

rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); light = data_ss(:,lightCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);
chosen_rref = 2.8; %quantile(rref,[0.5]);
deltar_grid = 0.5:0.25:3;%quantile(rref,[0.05 0.95])
y_grid = 0.075:0.025:0.35;%quantile(y,[0.05 0.95])
rdot_response1 = nan(length(y_grid),length(deltar_grid));
Amean_response2 = nan(length(y_grid),length(deltar_grid)); % from rdot model
% Amean_predicted = nan(length(y_grid),length(deltar_grid)); % from Amean model
% coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
% coeffs_Amean = [1.71368183  0.46295979  0.04482469  0.14650472  0.44938972 -1.53883634 -0.55715188 0];


for ct_y=1:length(y_grid)
    for ct_deltar=1:length(deltar_grid)
        y0 = y_grid(ct_y); t0 = 0; 
        thisDeltar = deltar_grid(ct_deltar);
        thisRref = chosen_rref;
        r0 = thisRref - thisDeltar;
        if r0 > 0
            V0 = r0*y0;
            tspan = [0:0.01:1.5];
            ystart = [y0 V0]; 
            % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
            log_rdot = coeffs(1) + coeffs(2)*log(y0) + coeffs(4) + coeffs(5)*log(thisDeltar) + ...
                coeffs(6)*log(thisRref) + coeffs(7)*log(y0)*log(thisDeltar) + ...
                + coeffs(8)*log(thisRref)*log(thisDeltar);
            rdot = exp(log_rdot);
            [t,y] = ode45(@(t,y) filteredState_BlindLandingtrack.const_drdt_movement(t,y,rdot), tspan, ystart);
            r = y(:,2)./y(:,1);
            tv_interp = interp1(r,[t y],thisRref);
            if tv_interp(2) > 0 && r(end) >= thisRref
                deltat = tv_interp(1); deltaV = tv_interp(3)-V0;
                Amean_response2(ct_y, ct_deltar) = deltaV/deltat;
                
                rdot_response1(ct_y, ct_deltar) = rdot;
                
%                 log_Amean_pred = coeffs_Amean(1) + coeffs_Amean(2)*log(y0) + coeffs_Amean(4) + coeffs_Amean(5)*log(thisDeltar) + ...
%                 coeffs_Amean(6)*log(chosen_rref) + coeffs_Amean(7)*log(y0)*log(thisDeltar) + ...
%                 + coeffs_Amean(8)*log(chosen_rref)*log(thisDeltar);
%                 Amean_predicted(ct_y, ct_deltar) = exp(log_Amean_pred);
            end
            
            
            
            
        else
            continue;
        end
    end
end

% Collect data points for plotting
binwidth = 0.3; 
data_ss2 = data_ss(data_ss(:,lightCol) == 3 & (data_ss(:,rrefCol) >= chosen_rref-binwidth/2 & data_ss(:,rrefCol) <= chosen_rref+binwidth/2),:);
rdot = data_ss2(:,rdotCol); y = data_ss2(:,yCol); light = data_ss2(:,lightCol); deltar = data_ss2(:,deltarCol); rref = data_ss2(:,rrefCol); amean = data_ss2(:,ameanCol);

% managing colorbar
% rdot_range = [8 25];
rdot_range = [5 20];
nColors = round(diff(rdot_range)/1);
rdot_cmap = jet(2*nColors);
rdot_cmap = rdot_cmap(round(nColors/2):round(nColors/2)+nColors-1,:);
rdot_edges = round(linspace(rdot_range(1),rdot_range(2),size(rdot_cmap,1)+1),2); % # edges = # color bins + 1

% amean_range = [0 10];
amean_range = [0 2.5];
nColors = round(diff(amean_range)/0.5);
amean_cmap = jet(2*nColors);
amean_cmap = amean_cmap(round(nColors/2):round(nColors/2)+nColors-1,:);
% amean_cmap = jet(round(diff(amean_range)/0.5));  
amean_edges = round(linspace(amean_range(1),amean_range(2),size(amean_cmap,1)+1),2); % # edges = # color bins + 1

[Y, DELTAR] = meshgrid(deltar_grid, y_grid);
rdot_fig = figure; hold on;
% subplot(1,2,1);
colormap(rdot_cmap);
[~,c1] = contourf(DELTAR, Y, rdot_response1, rdot_edges, 'LineStyle', 'none');

[data_cmap, ~] = discretize(rdot, rdot_edges);
data_cmap(rdot<=rdot_edges(1)) = 1; data_cmap(rdot>=rdot_edges(end)) = size(rdot_cmap,1);
scatter(y, deltar, 25, rdot_cmap(data_cmap',:),'filled','o','MarkerEdgeColor',[0 .5 .5]);

title('rdot (s-2)', 'FontSize', 14);
xlabel('y (m)', 'FontSize', 14);
ylabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);
% xlim([0.1 0.35]);
% xticks([0.1:0.05:0.35]);
% ylim([2 5]);
% yticks([2:0.5:5]);

cmap_bar = colorbar('eastoutside');
caxis(rdot_edges([1 end]));
cmap_bar.Ticks = [rdot_edges(1) rdot_edges(end)];
% cmap_bar.Label.String = 'rdot';
cmap_bar.Label.FontSize = 16;

amean_fig = figure; hold on;
% subplot(1,2,2) % rdot=f(y, deltar)
colormap(amean_cmap);
contourf(DELTAR, Y, Amean_response2, amean_edges, 'LineStyle', 'none');

[data_cmap, ~] = discretize(amean, amean_edges);
data_cmap(amean<=amean_edges(1)) = 1; data_cmap(amean>=amean_edges(end)) = size(amean_cmap,1);
scatter(y, deltar, 25, amean_cmap(data_cmap',:),'filled','o','MarkerEdgeColor',[0 .5 .5]);

title('Amean (ms-2)', 'FontSize', 14);
xlabel('y (m)', 'FontSize', 14);
ylabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);
% xlim([0.1 0.35]);
% xticks([0.1:0.05:0.35]);
% ylim([2 5]);
% yticks([2:0.5:5]);

cmap_bar = colorbar('eastoutside');
caxis(amean_edges([1 end]));
cmap_bar.Ticks = [amean_edges(1) amean_edges(end)];
% cmap_bar.Label.String = 'rdot';
cmap_bar.Label.FontSize = 16;

%% Mean accleration plot (TO BE USED) - with deltar and rref on two axes and for a fixed y
% close all;
% ameanCol = 17;
ameanCol = 18;
chosen_light = 3;

rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); light = data_ss(:,lightCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);
chosen_y = 0.22; %quantile(y,[0.5]);
deltar_grid = 0.5:0.5:7;%quantile(rref,[0.05 0.95])
rref_grid = 1:0.5:8;%quantile(y,[0.05 0.95])
rdot_response1 = nan(length(y_grid),length(deltar_grid));
Amean_response2 = nan(length(y_grid),length(deltar_grid)); % from rdot model
Amean_predicted = nan(length(y_grid),length(deltar_grid)); % from Amean model
% coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
coeffs_Amean = [1.58702290  0.42800399  0.05374526  0.16907431  0.32467186 -1.45318255 -0.58756146 0];


for ct_rref=1:length(rref_grid)
    for ct_deltar=1:length(deltar_grid)
        y0 = chosen_y; t0 = 0; 
        thisDeltar = deltar_grid(ct_deltar);
        thisRref = rref_grid(ct_rref);
        r0 = thisRref - thisDeltar;
        if r0 > 0
            V0 = r0*y0;
            tspan = [0:0.01:1.5];
            ystart = [y0 V0]; 
            % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
            log_rdot = coeffs(1) + coeffs(2)*log(y0) + coeffs(5)*log(thisDeltar) + ...
                coeffs(6)*log(thisRref) + coeffs(7)*log(y0)*log(thisDeltar) + ...
                + coeffs(8)*log(thisRref)*log(thisDeltar);
            if chosen_light == 3
                log_rdot = log_rdot + coeffs(4);
            end
            rdot = exp(log_rdot);
            [t,y] = ode45(@(t,y) filteredState_BlindLandingtrack.const_drdt_movement(t,y,rdot), tspan, ystart);
            r = y(:,2)./y(:,1);
            firstIndex = find(r>thisRref, 1);
            if isempty(firstIndex)
                keyboard;
            else
%                 r = r(1:firstIndex);
                tv_interp = interp1(r(1:firstIndex),[t(1:firstIndex) y(1:firstIndex, :)], thisRref);
%                 figHandle = figure; hold on;
%                 plot(t(1:firstIndex),r(1:firstIndex));
%                 yline(thisRref, 'r');
%                 yline(r0, 'g');
%                 ylim([0 thisRref+1]);
%                 pause(2);
%                 close(figHandle);
                
%             if r(end) <= thisRref
%                 keyboard;
            end
            if tv_interp(2) > 0 %&& r(end) >= thisRref
                deltat = tv_interp(1); deltaV = tv_interp(3)-V0;
%                 Amean_response2(ct_rref, ct_deltar) = deltaV/deltat;
                y = [y(1:firstIndex-1, :); tv_interp(2:3)];
                Amean_response2(ct_rref, ct_deltar) = mean((rdot*y(:,1).^2-y(:,2).^2)./y(:,1));
                
                rdot_response1(ct_rref, ct_deltar) = rdot;
                
                log_Amean_pred = coeffs_Amean(1) + coeffs_Amean(2)*log(y0) + coeffs_Amean(4) + coeffs_Amean(5)*log(thisDeltar) + ...
                coeffs_Amean(6)*log(thisRref) + coeffs_Amean(7)*log(y0)*log(thisDeltar) + ...
                + coeffs_Amean(8)*log(thisRref)*log(thisDeltar);
                Amean_predicted(ct_rref, ct_deltar) = exp(log_Amean_pred);
            end
            
            
            
            
        else
            continue;
        end
    end
end

% Collect data points for plotting
binwidth = 0.05; 
data_ss2 = data_ss(data_ss(:,lightCol) == chosen_light & (data_ss(:,yCol) >= chosen_y-binwidth/2 & data_ss(:,yCol) <= chosen_y+binwidth/2),:);
rdot = data_ss2(:,rdotCol); y = data_ss2(:,yCol); light = data_ss2(:,lightCol); deltar = data_ss2(:,deltarCol); rref = data_ss2(:,rrefCol); amean = data_ss2(:,ameanCol);

% managing colorbar
% rdot_range = [8 25];
rdot_range = [5 20];
nColors = round(diff(rdot_range)/1);
rdot_cmap = jet(2*nColors);
rdot_cmap = rdot_cmap(round(nColors/2):round(nColors/2)+nColors-1,:);
rdot_edges = round(linspace(rdot_range(1),rdot_range(2),size(rdot_cmap,1)+1),2); % # edges = # color bins + 1

% amean_range = [0 10];
amean_range = [0 2.5];
nColors = round(diff(amean_range)/0.5);
amean_cmap = jet(2*nColors);
amean_cmap = amean_cmap(round(nColors/2):round(nColors/2)+nColors-1,:);
% amean_cmap = jet(round(diff(amean_range)/0.5));  
amean_edges = round(linspace(amean_range(1),amean_range(2),size(amean_cmap,1)+1),2); % # edges = # color bins + 1

[RREF, DELTAR] = meshgrid(deltar_grid, rref_grid);
rdot_fig = figure; hold on;
% subplot(1,2,1);
colormap(rdot_cmap);
[~,c1] = contourf(DELTAR, RREF, rdot_response1, rdot_edges, 'LineStyle', 'none');

[data_cmap, ~] = discretize(rdot, rdot_edges);
data_cmap(rdot<=rdot_edges(1)) = 1; data_cmap(rdot>=rdot_edges(end)) = size(rdot_cmap,1);
scatter(rref, deltar, 25, rdot_cmap(data_cmap',:),'filled','o','MarkerEdgeColor',[0 .5 .5]);

title('rdot (s-2)', 'FontSize', 14);
xlabel('rref (s-1)', 'FontSize', 14);
ylabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);
% xlim([0.1 0.35]);
% xticks([0.1:0.05:0.35]);
% ylim([2 5]);
% yticks([2:0.5:5]);

cmap_bar = colorbar('eastoutside');
caxis(rdot_edges([1 end]));
cmap_bar.Ticks = [rdot_edges(1) rdot_edges(end)];
% cmap_bar.Label.String = 'rdot';
cmap_bar.Label.FontSize = 16;

amean_fig = figure; hold on;
% subplot(1,2,2) % rdot=f(y, deltar)
colormap(amean_cmap);
contourf(DELTAR, RREF, Amean_response2, amean_edges, 'LineStyle', 'none');
% contourf(DELTAR, RREF, Amean_predicted, amean_edges, 'LineStyle', 'none');

[data_cmap, ~] = discretize(amean, amean_edges);
data_cmap(amean<=amean_edges(1)) = 1; data_cmap(amean>=amean_edges(end)) = size(amean_cmap,1);
scatter(rref, deltar, 25, amean_cmap(data_cmap',:),'filled','o','MarkerEdgeColor',[0 .5 .5]);

title('Amean (ms-2)', 'FontSize', 14);
xlabel('rref (s-1)', 'FontSize', 14);
ylabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);
% xlim([0.1 0.35]);
% xticks([0.1:0.05:0.35]);
% ylim([2 5]);
% yticks([2:0.5:5]);

cmap_bar = colorbar('eastoutside');
caxis(amean_edges([1 end]));
cmap_bar.Ticks = [amean_edges(1) amean_edges(end)];
% cmap_bar.Label.String = 'rdot';
cmap_bar.Label.FontSize = 16;

figure;
histogram([Amean_response2(:)])
xlabel('Simulated Amean (ms-2)', 'FontSize', 14);
figure;
histogram(amean);
xlabel('Observed Amean (ms-2)', 'FontSize', 14);

% Comparison of rdot model and Amean model
figure;
subplot(1,2,1) 
contourf(DELTAR, RREF, Amean_predicted, 10,'ShowText','on');
title('Amean from Amean model (ms-2) (s-2)', 'FontSize', 14);
xlabel('rref (s-1)', 'FontSize', 14);
ylabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);

subplot(1,2,2) % rdot=f(y, deltar)
contourf(DELTAR, RREF, Amean_response2, 10,'ShowText','on');
% contourf(Amean_response2, Y, DELTAR, 10,'ShowText','on');
title('Amean from rdot model (ms-2)', 'FontSize', 14);
xlabel('rref (s-1)', 'FontSize', 14);
ylabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);

%% Comparison of Amean in different light conditions - fixed y and varying rref and deltar
ameanCol = 18;

rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); light = data_ss(:,lightCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);
chosen_y = 0.22; %quantile(y,[0.5]);
deltar_grid = 0.5:0.5:7;%quantile(rref,[0.05 0.95])
rref_grid = 1:0.5:8;%quantile(y,[0.05 0.95])
light_grid = 1:1:3;
rdot_response1 = nan(length(y_grid),length(deltar_grid), length(light_grid));
Amean_response2 = nan(length(y_grid),length(deltar_grid), length(light_grid)); % from rdot model
% Amean_predicted = nan(length(y_grid),length(deltar_grid)); % from Amean model
% coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
% coeffs_Amean = [1.71368183  0.46295979  0.04482469  0.14650472  0.44938972 -1.53883634 -0.55715188 0];

for ct_light = 1:length(light_grid)
    for ct_rref=1:length(rref_grid)
        for ct_deltar=1:length(deltar_grid)
            thisLight = light_grid(ct_light);
            y0 = chosen_y; t0 = 0; 
            thisDeltar = deltar_grid(ct_deltar);
            thisRref = rref_grid(ct_rref);
            r0 = thisRref - thisDeltar;
            if r0 > 0
                V0 = r0*y0;
                tspan = [0:0.01:1.5];
                ystart = [y0 V0]; 
                % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                log_rdot = coeffs(1) + coeffs(2)*log(y0) + coeffs(5)*log(thisDeltar) + ...
                    coeffs(6)*log(thisRref) + coeffs(7)*log(y0)*log(thisDeltar) + ...
                    + coeffs(8)*log(thisRref)*log(thisDeltar);
                if thisLight == 2
                    log_rdot = log_rdot + coeffs(3);
                elseif thisLight == 3
                    log_rdot = log_rdot + coeffs(4);
                end
                rdot = exp(log_rdot);
                [t,y] = ode45(@(t,y) filteredState_BlindLandingtrack.const_drdt_movement(t,y,rdot), tspan, ystart);
                r = y(:,2)./y(:,1);
                firstIndex = find(r>thisRref, 1);
                if isempty(firstIndex)
                    keyboard;
                else
    %                 r = r(1:firstIndex);
                    tv_interp = interp1(r(1:firstIndex),[t(1:firstIndex) y(1:firstIndex, :)], thisRref);
    %             if r(end) <= thisRref
    %                 keyboard;
                end
                if tv_interp(2) > 0 %&& r(end) >= thisRref
                    deltat = tv_interp(1); deltaV = tv_interp(3)-V0;
                    Amean_response2(ct_rref, ct_deltar, ct_light) = deltaV/deltat;

                    rdot_response1(ct_rref, ct_deltar, ct_light) = rdot;

    %                 log_Amean_pred = coeffs_Amean(1) + coeffs_Amean(2)*log(y0) + coeffs_Amean(4) + coeffs_Amean(5)*log(thisDeltar) + ...
    %                 coeffs_Amean(6)*log(chosen_rref) + coeffs_Amean(7)*log(y0)*log(thisDeltar) + ...
    %                 + coeffs_Amean(8)*log(chosen_rref)*log(thisDeltar);
    %                 Amean_predicted(ct_y, ct_deltar) = exp(log_Amean_pred);
                end




            else
                continue;
            end
        end
    end
end

%% Comparison of Amean in different light conditions - median values of y, rref and deltar
rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); light = data_ss(:,lightCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);
chosen_y = quantile(y,[0.5]);
chosen_deltar = quantile(deltar,[0.5]);
chosen_rref = quantile(rref,[0.5]);
light_grid = 1:1:3;
rdot_response1 = nan(length(light_grid),1);
Amean_response2 = nan(length(light_grid),1); % from rdot model

for ct_light = 1:length(light_grid)
    
    thisLight = light_grid(ct_light);
    y0 = chosen_y; t0 = 0;
    thisDeltar = chosen_deltar;
    thisRref = chosen_rref;
    r0 = thisRref - thisDeltar;
    if r0 > 0
        V0 = r0*y0;
        tspan = [0:0.01:1.5];
        ystart = [y0 V0];
        % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
        log_rdot = coeffs(1) + coeffs(2)*log(y0) + coeffs(5)*log(thisDeltar) + ...
            coeffs(6)*log(thisRref) + coeffs(7)*log(y0)*log(thisDeltar) + ...
            + coeffs(8)*log(thisRref)*log(thisDeltar);
        if thisLight == 2
            log_rdot = log_rdot + coeffs(3);
        elseif thisLight == 3
            log_rdot = log_rdot + coeffs(4);
        end
        rdot = exp(log_rdot);
        [t,y] = ode45(@(t,y) filteredState_BlindLandingtrack.const_drdt_movement(t,y,rdot), tspan, ystart);
        r = y(:,2)./y(:,1);
        firstIndex = find(r>thisRref, 1);
        if isempty(firstIndex)
            keyboard;
        else
            %                 r = r(1:firstIndex);
            tv_interp = interp1(r(1:firstIndex),[t(1:firstIndex) y(1:firstIndex, :)], thisRref);
            %             if r(end) <= thisRref
            %                 keyboard;
        end
        if tv_interp(2) > 0 %&& r(end) >= thisRref
            deltat = tv_interp(1); deltaV = tv_interp(3)-V0;
            Amean_response2(ct_light) = deltaV/deltat;
            
            rdot_response1(ct_light) = rdot;
            
            %                 log_Amean_pred = coeffs_Amean(1) + coeffs_Amean(2)*log(y0) + coeffs_Amean(4) + coeffs_Amean(5)*log(thisDeltar) + ...
            %                 coeffs_Amean(6)*log(chosen_rref) + coeffs_Amean(7)*log(y0)*log(thisDeltar) + ...
            %                 + coeffs_Amean(8)*log(chosen_rref)*log(thisDeltar);
            %                 Amean_predicted(ct_y, ct_deltar) = exp(log_Amean_pred);
        end



        
    else
        continue;
    end
end

amean_light1 = data_ss(data_ss(:,lightCol) == 1,ameanCol);
amean_light2 = data_ss(data_ss(:,lightCol) == 2,ameanCol);
amean_light3 = data_ss(data_ss(:,lightCol) == 3,ameanCol);

map = brewermap(3,'Set1'); 
figure;
histogram(amean_light1,'facecolor',map(1,:),'facealpha',.5,'edgecolor','none');
hold on;
% xline(median(amean_light1), 'LineColor',map(1,:));

% histogram(amean_light3,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none');
% hold on;
histogram(amean_light3,'facecolor',map(3,:),'facealpha',.5,'edgecolor','none');
legend('low','high','fontsize',16)
xlabel('Amean (ms-2)', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);
axis tight

figure;
subplot(1,2,1);
histfit(amean_light1,[],'Gamma')
pd = fitdist(amean_light1,'Gamma');
subplot(1,2,2);
histfit(amean_light3,[],'Gamma')
pd = fitdist(amean_light3,'Gamma');

[median(amean_light1) Amean_response2(1)]
[median(amean_light2) Amean_response2(2)]
[median(amean_light3) Amean_response2(3)]
% Collect data points for plotting
binwidth = 0.4; 
data_ss2 = data_ss(data_ss(:,lightCol) == chosen_light & (data_ss(:,yCol) >= chosen_y-binwidth/2 & data_ss(:,yCol) <= chosen_y+binwidth/2),:);
rdot = data_ss2(:,rdotCol); y = data_ss2(:,yCol); light = data_ss2(:,lightCol); deltar = data_ss2(:,deltarCol); rref = data_ss2(:,rrefCol); amean = data_ss2(:,ameanCol);

% managing colorbar
% rdot_range = [8 25];
rdot_range = [5 20];
nColors = round(diff(rdot_range)/1);
rdot_cmap = jet(2*nColors);
rdot_cmap = rdot_cmap(round(nColors/2):round(nColors/2)+nColors-1,:);
rdot_edges = round(linspace(rdot_range(1),rdot_range(2),size(rdot_cmap,1)+1),2); % # edges = # color bins + 1

% amean_range = [0 10];
amean_range = [0 2.5];
nColors = round(diff(amean_range)/0.5);
amean_cmap = jet(2*nColors);
amean_cmap = amean_cmap(round(nColors/2):round(nColors/2)+nColors-1,:);
% amean_cmap = jet(round(diff(amean_range)/0.5));  
amean_edges = round(linspace(amean_range(1),amean_range(2),size(amean_cmap,1)+1),2); % # edges = # color bins + 1

[RREF, DELTAR] = meshgrid(deltar_grid, rref_grid);
rdot_fig = figure; hold on;
% subplot(1,2,1);
colormap(rdot_cmap);
[~,c1] = contourf(DELTAR, RREF, rdot_response1, rdot_edges, 'LineStyle', 'none');

[data_cmap, ~] = discretize(rdot, rdot_edges);
data_cmap(rdot<=rdot_edges(1)) = 1; data_cmap(rdot>=rdot_edges(end)) = size(rdot_cmap,1);
scatter(rref, deltar, 25, rdot_cmap(data_cmap',:),'filled','o','MarkerEdgeColor',[0 .5 .5]);

title('rdot (s-2)', 'FontSize', 14);
xlabel('y (m)', 'FontSize', 14);
ylabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);
% xlim([0.1 0.35]);
% xticks([0.1:0.05:0.35]);
% ylim([2 5]);
% yticks([2:0.5:5]);

cmap_bar = colorbar('eastoutside');
caxis(rdot_edges([1 end]));
cmap_bar.Ticks = [rdot_edges(1) rdot_edges(end)];
% cmap_bar.Label.String = 'rdot';
cmap_bar.Label.FontSize = 16;

amean_fig = figure; hold on;
% subplot(1,2,2) % rdot=f(y, deltar)
colormap(amean_cmap);
contourf(DELTAR, RREF, Amean_response2, amean_edges, 'LineStyle', 'none');

[data_cmap, ~] = discretize(amean, amean_edges);
data_cmap(amean<=amean_edges(1)) = 1; data_cmap(amean>=amean_edges(end)) = size(amean_cmap,1);
scatter(rref, deltar, 25, amean_cmap(data_cmap',:),'filled','o','MarkerEdgeColor',[0 .5 .5]);

title('Amean (ms-2)', 'FontSize', 14);
xlabel('rref (s-1)', 'FontSize', 14);
ylabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);
% xlim([0.1 0.35]);
% xticks([0.1:0.05:0.35]);
% ylim([2 5]);
% yticks([2:0.5:5]);

cmap_bar = colorbar('eastoutside');
caxis(amean_edges([1 end]));
cmap_bar.Ticks = [amean_edges(1) amean_edges(end)];
% cmap_bar.Label.String = 'rdot';
cmap_bar.Label.FontSize = 16;


%% Mean accleration plot (NOT TO BE USED) - with deltar and y on two axes and for a fixed rref
rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); light = data_ss(:,lightCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);
chosen_rref = quantile(rref,[0.75]);
deltar_grid = 0.5:0.05:3.7;%quantile(deltar,[0.05 0.95])
y_grid = 0.05:0.025:0.35;%quantile(y,[0.05 0.95])
rdot_response1 = nan(length(y_grid),length(deltar_grid));
Amean_response2 = nan(length(y_grid),length(deltar_grid)); % from rdot model
Amean_predicted = nan(length(y_grid),length(deltar_grid)); % from Amean model
% coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
coeffs_Amean = [1.58702290  0.42800399  0.05374526  0.16907431  0.32467186 -1.45318255 -0.58756146 0];


for ct_y=1:length(y_grid)
    for ct_deltar=1:length(deltar_grid)
        y0 = y_grid(ct_y); t0 = 0; thisDeltar = deltar_grid(ct_deltar);
        r0 = chosen_rref - thisDeltar;
        if r0 > 0
            V0 = r0*y0;
            tspan = [0:0.01:1];
            ystart = [y0 V0]; 
            % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
            log_rdot = coeffs(1) + coeffs(2)*log(y0) + coeffs(4) + coeffs(5)*log(thisDeltar) + ...
                coeffs(6)*log(chosen_rref) + coeffs(7)*log(y0)*log(thisDeltar) + ...
                + coeffs(8)*log(chosen_rref)*log(thisDeltar);
            rdot = exp(log_rdot);
            [t,y] = ode45(@(t,y) filteredState_BlindLandingtrack.const_drdt_movement(t,y,rdot), tspan, ystart);
            r = y(:,2)./y(:,1);
            tv_interp = interp1(r,[t y],chosen_rref);
            if tv_interp(2) > 0 && r(end) >= chosen_rref
                deltat = tv_interp(1); deltaV = tv_interp(3)-V0;
                Amean_response2(ct_y, ct_deltar) = deltaV/deltat;
                
                log_Amean_pred = coeffs_Amean(1) + coeffs_Amean(2)*log(y0) + coeffs_Amean(4) + coeffs_Amean(5)*log(thisDeltar) + ...
                coeffs_Amean(6)*log(chosen_rref) + coeffs_Amean(7)*log(y0)*log(thisDeltar) + ...
                + coeffs_Amean(8)*log(chosen_rref)*log(thisDeltar);
                Amean_predicted(ct_y, ct_deltar) = exp(log_Amean_pred);
            end
            
            rdot_response1(ct_y, ct_deltar) = rdot;
            
            
        else
            continue;
        end
    end
end
[Y, DELTAR] = meshgrid(deltar_grid, y_grid);
figure;
subplot(1,2,1) 
contourf(DELTAR, Y, rdot_response1, 10,'ShowText','on');
title('rdot (s-2)', 'FontSize', 14);
xlabel('y (m)', 'FontSize', 14);
ylabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);

subplot(1,2,2) % rdot=f(y, deltar)
contourf(DELTAR, Y, Amean_response2, 10,'ShowText','on');
% contourf(Amean_response2, Y, DELTAR, 10,'ShowText','on');
title('Amean (ms-2)', 'FontSize', 14);
xlabel('y (m)', 'FontSize', 14);
ylabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);


figure;
subplot(1,2,1) 
contourf(DELTAR, Y, Amean_response2, 10,'ShowText','on');
% contourf(Amean_response2, Y, DELTAR, 10,'ShowText','on');
title('Amean (ms-2) from rdot model', 'FontSize', 14);
xlabel('y (m)', 'FontSize', 14);
ylabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);

subplot(1,2,2) % rdot=f(y, deltar)
contourf(DELTAR, Y, Amean_predicted, 10,'ShowText','on');
% contourf(Amean_response2, Y, DELTAR, 10,'ShowText','on');
title('Amean (ms-2) from Amean model', 'FontSize', 14);
xlabel('y (m)', 'FontSize', 14);
ylabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);

%% Plot rdot as a function of y, delta_r and rref 
% y as x1, deltar as x2 and rref as x3
close all;
rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); light = data_ss(:,lightCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);

% Untransformed domain
res = rdot; x1 = y; x2 = deltar; x3 = rref; x4 = light;

% plot3(log(rref), log(rdot0), log(delta_t),'.')

x3_centers = quantile(x3,[0.25 0.5 0.75]); % vlaues for which rref data is binned
binwidth = 1.5; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)

x2_forplot = [];
for ct=1:3 % for each light condition
    for ct1=1:length(x3_centers)  % for each x3
%         data_ss2 = data_ss(data_ss(:,lightCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2),:);
%         deltar_ss2 = data_ss2(:,deltarCol);
%         [min(deltar_ss2) max(deltar_ss2)]
%         
        dummy = x2(data_ss(:,lightCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2));
%         dummy1 = data_ss(data_ss(:,lightCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2), deltarCol);
%         [min(dummy) max(dummy)]
        x2_forplot = [x2_forplot; dummy]; 
    end
end
x2range = [min(x2_forplot) max(x2_forplot)];
x2range = [round(x2range(1)*10)/10 round(max(x2range(2))*10)/10];

% cmap=brewermap(round(diff(x2range)/0.5),'Set1'); 
% cmap = jet(round(diff(x2range)/0.25));  
cmap = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880];
edges = round(linspace(x2range(1),x2range(2),size(cmap,1)+1),2); % # edges = # color bins + 1

% For plotting continuous lines
x2quantiles = quantile(x2_forplot,[0.05 0.5 0.95]) % deltar values for which lines are plotted
% x2quantiles = quantile(data_ss(:,deltarCol),[0.05 0.5 0.95]) % deltar values for which lines are plotted
% x2quantiles = [0.6192 1.7449 3.79];
y_values_cont = 0.05:0.01:0.35;

figure2d = figure; hold on;
colormap(cmap);

figure2d_log = figure; hold on;
colormap(cmap);
for ct=1:3 % for each light condition
    for ct1=1:length(x3_centers)  % for each x3
        figure(figure2d);
        subplot(length(x3_centers), 3,(ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,lightCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2),:);
        deltar_ss2 = data_ss2(:,deltarCol);
%         [min(deltar_ss2) max(deltar_ss2)]
        [data_cmap, ~] = discretize(deltar_ss2, edges);
        data_cmap(deltar_ss2<=edges(1)) = 1; data_cmap(deltar_ss2>=edges(end)) = size(cmap,1);
        scatter(data_ss2(:,yCol), data_ss2(:,rdotCol), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
%         x2quantiles = quantile(deltar_ss2,[0.05 0.5 0.95]) % deltar values for which lines are plotted
%         dum = [];
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(y_values_cont) + coeffs(5)*log(x2quantiles(ct2)) + ...
                coeffs(6)*log(x3_centers(ct1)) + coeffs(7)*log(y_values_cont)*log(x2quantiles(ct2)) + ...
                + coeffs(8)*log(x3_centers(ct1))*log(x2quantiles(ct2));
            if ct == 2
                % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif ct == 3
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
%             dum(ct2) = linecolor;
            plot(y_values_cont, exp(rdot_response),'Color',cmap(linecolor,:),'LineWidth', 2);
        end
%         display(dum)
        
        set(gca, 'FontSize', 16);
        ylim([0, 40]);
        
        figure(figure2d_log);
        subplot(length(x3_centers), 3,(ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,lightCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2),:);
        deltar_ss2 = data_ss2(:,deltarCol);
%         [min(deltar_ss2) max(deltar_ss2)]
        [data_cmap, ~] = discretize(deltar_ss2, edges);
        data_cmap(deltar_ss2<=edges(1)) = 1; data_cmap(deltar_ss2>=edges(end)) = size(cmap,1);
        scatter(log(data_ss2(:,yCol)), log(data_ss2(:,rdotCol)), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
        x2quantiles = quantile(deltar_ss2,[0.05 0.5 0.95]); % deltar values for which lines are plotted
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(y_values_cont) + coeffs(5)*log(x2quantiles(ct2)) + ...
                coeffs(6)*log(x3_centers(ct1)) + coeffs(7)*log(y_values_cont)*log(x2quantiles(ct2)) + ...
                + coeffs(8)*log(x3_centers(ct1))*log(x2quantiles(ct2));
            if ct == 2
                % coeffs = [intercept logy light2 light3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif ct == 3
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            plot(log(y_values_cont), rdot_response,'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        
        set(gca, 'FontSize', 16);
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
cmap_bar.Label.String = 'delta_r';
cmap_bar.Label.FontSize = 16;
ylabel('rdot (s-2)', 'FontSize', 14);
xlabel('y (m)', 'FontSize', 14);
%% Functions used in this script
function figHandle = createBoxPlot(variable, labels, yxislabel)
    % for boxplots of a cell array
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

    
    figHandle = figure; hold on;
    figure(figHandle);
    boxplot2(variable, 'Labels', labels, 'OutlierSize', 0.00001);
%     hbox = gca;
    set(gca, 'FontSize', 18); %grid on;
%     xlabel('Treatments', 'FontSize', 18);
    ylabel(yxislabel, 'FontSize', 18);
    % ylim([0 1.5]);
    for ct = 1:length(variable)
            x=(ct+(rand(length(variable{ct, 1}),1)-0.5)/4);

            f = scatter(x(:,1),variable{ct, 1},10,'k','filled'); 
            f.MarkerFaceAlpha = 0.5;
    %         keyboard;
    end
end






