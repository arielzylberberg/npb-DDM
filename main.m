addpath('./functions/');
addpath(genpath('./bads-master'));

%%

load('./fake_data/fake_data_exp');

d = true_params;

%%
pars.error_flag = 2; %1 is perf, 2 is likelihood
pars.bandwidth = 0.05;
plot_flag = 1;

%% just infer the bounds, assuming the other parameters are known

theta(1) = d.kappa;
theta(2) = d.ndt_mu;
theta(3) = d.ndt_sigma;
theta(4) = 0;
theta(5) = 0;


[~,P] = wrapper_dtb_empiricalbound_rt(theta,rt,coh,choice,c,pars,plot_flag);

figure
plot(d.t,d.Bup,P.t,P.Bup);
hold all
plot(d.t,-1*d.Bup,P.t,-1*P.Bup);
xlabel('Time [s]');
ylabel('Accumulated evidence [a.u.]');
legend('Ground truth','recovered');

%% fitting

% kappa, ndt_mu, ndt_sigma, coh0, y0
tl = [5,  0.1, .001 ,0,0];
th = [40, 0.7, .15 ,0,0];
tg = [15, 0.2, .02 ,0,0];


fn_fit = @(theta) (wrapper_dtb_empiricalbound_rt(theta,rt,coh,choice,c,pars,plot_flag));

% fit
options = optimset('Display','final','TolFun',.1,'FunValCheck','on');
ptl = tl;
pth = th;
[theta,fval,exitflag,output] = bads(@(theta) fn_fit(theta),tg,tl,th,ptl,pth,options);

% eval best
[err,P] = wrapper_dtb_empiricalbound_rt(theta,rt,coh,choice,c,pars,plot_flag);

% plot inferred bounds against ground truth
figure
plot(d.t,d.Bup,P.t,P.Bup);
hold all
plot(d.t,-1*d.Bup,P.t,-1*P.Bup);
xlabel('Time [s]');
ylabel('Accumulated evidence [a.u.]');
legend('Ground truth','recovered');






