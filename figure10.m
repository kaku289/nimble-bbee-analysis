%% Figure10
% Plotting data from https://royalsocietypublishing.org/doi/pdf/10.1098/rspb.2020.3051
Vg = [379 356 362 345 300]./1000; % Mean ground speed in different winds for cross-hatch pattern
h = [87 93 113 134 142]./1000;
winds = [2 1 0 -1 -2]; % positive winds are difficult ones
setpoints = Vg./h;
airspeed = Vg-winds;
hbee_fit = fitlm(airspeed,setpoints);

%% bbee data (my own sweet wonderful amazing data)
rref = [2.07, 2.08, 2.09, 2.16, 2.24, 2.4]
winds = [0, 0.28, 0.98, 1.57, 2.54, 3.41];
bbee_fit = fitlm(winds,rref);
hold on;
plot(bbee_fit)