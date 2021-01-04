syms r(t) d w k rref r0 rdot0
Dr = diff(r)

ode = diff(r,t,2) + 2*d*w*Dr + w^2*r == k*w^2*rref
cond1 = r(0) == r0
cond2 = Dr(0) == rdot0
conds = [cond1 cond2]

rSol(t) = dsolve(ode,conds)
simplify(rSol(t))

rSol1(t) = subs(rSol(t),{k, rref, r0, rdot0},{1, 1, 0, 0})
simplify(rSol1(t))

(exp(-t*w*(d + (d^2 - 1)^(1/2)))*(d - (d^2 - 1)^(1/2)))/(2*(d^2 - 1)^(1/2)) - (exp(-t*w*(d - (d^2 - 1)^(1/2)))*(d + (d^2 - 1)^(1/2)))/(2*(d^2 - 1)^(1/2)) + 1


%%
syms t(r) d w k rref r0 rdot0
Dr = 1/diff(t)

ode = 1/diff(t,r,2) + 2*d*w*Dr + w^2*r == k*w^2*rref
cond1 = t(r0) == 0
cond2 = Dr(r0) == rdot0
conds = [cond1 cond2]

rSol(t) = dsolve(ode,conds)
simplify(rSol(t))

%%
tf = tfs; % load a transfer function. Here I selected track 151 from spokes, high condition
sysiddata % data used for identification

sys = ss(tf); % arbitrary x1 and x2
[ocsys, T] = canon(ss(tf),'companion'); % observable canonical form
ccsys = ss(ocsys.A', ocsys.C', ocsys.B', ocsys.D'); % contollable canonical form (https://nl.mathworks.com/help/ident/ug/canonical-state-space-realizations.html)
To2c = [ccsys.B ccsys.A*ccsys.B]/[ocsys.B ocsys.A*ocsys.B]; % state space coordinate transformation matrix from observable canonical form to controllable canonical form (http://controleducation.group.shef.ac.uk/statespace/state%20space%20feedback%203%20-%20transformation%20to%20get%20a%20canonical%20form.pdf)

[rsim,fit,x0] = compare(sysiddata{1}, tf)
x0_oc = T*x0
x0_cc = To2c*x0_oc

% For zero initial condition and different rref

rref = [1:0.5:8]';
t=0:0.005:0.5;
tr= nan*zeros(length(rref),1);
for ct_rref=1:length(rref)
    rsim = lsim(ccsys, rref(ct_rref)*ones(length(t),1), t, [0; 0]);
%     assert(rsim(end) > 0.9*rref(ct_rref))
    
    % find rise time
    indx = find(rsim>0.9*rref(ct_rref),1,'first');
    tr(ct_rref) = diff(interp1(rsim(1:indx+4),t(1:indx+4),[0.5 0.9]*rref(ct_rref)));
    
    
end

figure;
scatter(rref,tr,20,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
xlabel('r* (s-1)', 'FontSize', 18);
ylabel('t_r (s)', 'FontSize', 18);
          

% for non-zero r and zer rdot as initial condition and different rref
rref = [2.5:0.5:8]';
t=0:0.005:0.5;
tr= nan*zeros(length(rref),1);
for ct_rref=1:length(rref)
    rsim = lsim(ccsys, rref(ct_rref)*ones(length(t),1), t, [0.9; 25]);
%     assert(rsim(end) > 0.9*rref(ct_rref))
    
    % find rise time
    indx = find(rsim>0.9*rref(ct_rref),1,'first');
    tr(ct_rref) = diff(interp1(rsim(1:indx+4),t(1:indx+4),[0.5 0.9]*rref(ct_rref)));
end

figure;
scatter(log(rref),log(tr),20,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
xlabel('log(r*) (s-1)', 'FontSize', 18);
ylabel('log(t_r)(s)', 'FontSize', 18);          
[ones(length(rref),1) log(rref)]\log(tr)

