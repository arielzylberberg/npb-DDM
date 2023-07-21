function [err,P] = wrapper_dtb_empiricalbound_rt(theta,rt,coh,choice,c,pars,plot_flag)
% function [err,P] = wrapper_dtb_empiricalbound_rt(theta,rt,coh,choice,c,pars,plot_flag)

%  INPUT
    % theta(1): kappa (signal-to-noise)
    % theta(2): mean of non-decision time
    % theta(3): s.d. of non-decision time
    % theta(4): bias term, added to the trial's coherence
    % theta(5): bias term, as a starting point offset

% 2014: Ariel Zylberberg wrote it



error_flag = pars.error_flag;


if isfield(pars,'bandwidth')
    bandwidth = pars.bandwidth;
else
    bandwidth = 0.05;
end

% sanity check
if error_flag==3 && any(coh<0)
    error('no negative coherences for acc fit')
end

%%
kappa  = theta(1);
ndt_m  = theta(2);
ndt_s  = theta(3);
coh0   = theta(4);
y0a    = theta(5);

%%

if ~isempty(pars) && isfield(pars,'t')
    t = pars.t;
    dt = t(2)-t(1);
else
    dt = 0.0005;
    t  = 0:dt:10;
end

if ~isempty(pars) && isfield(pars,'y')
    y = pars.y;
else
    y  = linspace(-3,3,1500)';%esto puede ser problematico
%     y  = linspace(min(Blo)-0.3,max(Bup)+0.3,1500)';%esto puede ser problematico
end
y0 = zeros(size(y));
y0(findclose(y,0)) = 1;
y0 = y0/sum(y0);



%%

%ignore ndt_s for dect calculation
dect = rt - ndt_m;
inds = dect>0;
q = ksdensity(dect(inds),t,'kernel','epanechnikov','support','positive','bandwidth',bandwidth,'support','positive');


F = @dtb_fp_cc_searchbnd;


prior = Rtable(coh,1)/sum(Rtable(coh,1));


%%

drift = kappa * unique(coh + coh0);

notabs_flag = false;
P = feval(F,drift,t,prior,q,y,y0,notabs_flag);


if error_flag==1
    %% eval error, binomial
    I = [choice==1, choice==0];
    n = sum(I); %all trials
    y = [sum(c(I(:,1))),sum(c(I(:,2)))]; %correct
    aux = P.up.p.*prior;
    ucoh = unique(coh);
    yfit(1) = sum(aux.*(ucoh>0)+0.5*aux.*(ucoh==0))/sum(aux);
    aux = P.lo.p.*prior;
    yfit(2) = sum(aux.*(ucoh<0)+0.5*aux.*(ucoh==0))/sum(aux);
    yfit = yfit.*n;
    err = -2 * sum(log(binopdf(y,n,yfit./n)));
    
elseif error_flag==2
    
    %% likelihood

    [nlogl,pPred,upRT,loRT] = logl_choiceRT_1d(P,choice,rt,coh,ndt_m,ndt_s);


    err = nlogl / length(rt); % return as per trial
    
    
else
    err = nan;
    
end

%%
%% print
fprintf('err=%.3f kappa=%.2f ndt_mu=%.2f ndt_s=%.2f coh0=%.2f y0=%.2f \n',...
    err,kappa,ndt_m,ndt_s,coh0,y0a);

%%
if plot_flag
    m = prctile(rt,99.5);
    
    figure(1);clf
    
    subplot(3,1,1);
    plot(t,P.Bup,'k');
    hold all
    plot(t,P.Blo,'k');
    xlabel('Time [s]');
    ylabel('Accumulated evidence [a.u.]');
%     title('Inferred bounds');
    xlim([0,m]);
    
    subplot(3,1,2);
    if error_flag==3
        [tt,xx] = curva_media(c,coh,[],0);
    else
        [tt,xx] = curva_media(choice,coh,[],0);
    end
    plot(tt,xx,'b.-');
    hold all
    ucoh = nanunique(coh);
    plot(ucoh,P.up.p,'r-');
    legend('data','model');
    xlabel('Motion Coherence');
    ylabel('P. rightward');
    
    
    subplot(3,1,3);
    rt_model = (P.up.mean_t.*P.up.p+P.lo.mean_t.*P.lo.p)./(P.up.p+P.lo.p) + ndt_m; %is not exact because it
    % doesn't take into account the curtailing of the non-decision time
    % distribution
    [tt,xx,ss] = curva_media(rt,coh,[],0);
    terrorbar(tt,xx,ss,'marker','.','linestyle','-','color','b');
    hold all
    plot(ucoh,rt_model,'r-');
    legend('data','model');
    xlabel('Motion Coherence');
    ylabel('Response Time [s]');
    
    set(gcf,'Position',[576   135   394   730]);
    format_figure(gcf);
    
    drawnow
    
end