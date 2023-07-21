function P = dtb_fp_cc_searchbnd(drift,t,prior,dtdist,y,y0,notabs_flag)
% searches empirical bounds with chang cooper method
% - no hazard
% 2014: Ariel Zylberberg wrote it 

% sanity check
if length(prior)~=length(drift)
    error('error 1, wrong sizes')
end

%norm to 1
prior = prior(:)'/sum(prior);%column vec
defective_dt_dist = 1-cumsum(dtdist/nansum(dtdist));

dt = t(2) - t(1);
dy = y(2) - y(1);
ny = length(y);
nt = length(t);

Bup = nan(nt,1);%bound, assumed symmetric

P = struct('drift',drift,'t',t,'prior',prior,'dtdist',dtdist,'y',y,'y0',y0,...
    'notabs_flag',notabs_flag);

nd = length(drift);

% Preallocate
P.up.pdf_t = zeros(nd,nt);
P.lo.pdf_t = zeros(nd,nt);
if notabs_flag
    P.notabs.pdf = zeros(nd,ny,nt);
end
p_threshold = 1.0E-5; % Threshold for proportion un-terminated to stop simulation


M = chang_cooper_sparsematrix(drift,nd,ny,dy,dt);

yr = repmat(y(:),nd,1);
u = repmat(y0(:),nd,1);

for k = 2:nt
    
    u = M\u;
    
    ur = reshape(u,ny,nd);
    
    %search B(k)
    fun = @(x) (defective_dt_dist(k) - sum(sum(ur(y<x & y>-x,:),1).*prior));

    X0 = [0,max(y)];
    if sign(fun(X0(1)))~=sign(fun(X0(2)))
        [Bup(k),~,exitflag] = fzero(fun,X0);
    else
        Bup(k) = nan;
    end
        
    % Select density that has crossed bounds
    P.up.pdf_t(:,k) = sum(ur(y>=Bup(k),:),1);
    P.lo.pdf_t(:,k) = sum(ur(y<=-Bup(k),:),1);
    
    % Keep only density within bounds
    outofbounds = yr<=-Bup(k) | yr>=Bup(k);
    u(outofbounds) = 0;

    % Save if requested
    if notabs_flag
        P.notabs.pdf(:,:,k) = ur';
    end
    
    if sum(sum(ur,1)<p_threshold)==nd
        break;
    end
    
end

if notabs_flag
    P.notabs.pos_t = sum(P.notabs.pdf(:,y'>=0,:),2);
    P.notabs.neg_t = sum(P.notabs.pdf(:,y'< 0,:),2);
end

P.up.p = sum(P.up.pdf_t,2);
P.lo.p = sum(P.lo.pdf_t,2);

t = t(:);

P.up.mean_t = transpose(t'*P.up.pdf_t')./P.up.p;
P.lo.mean_t = transpose(t'*P.lo.pdf_t')./P.lo.p;

P.up.cdf_t = cumsum(P.up.pdf_t,2);
P.lo.cdf_t = cumsum(P.lo.pdf_t,2);

P.Bup = Bup;
P.Blo = -1*Bup;

