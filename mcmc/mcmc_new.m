function [models] = mcmc_new(Niter, stepsize, fun, x0, bounds,A, b, bcut, k, exit_cond,fail, varargin)
%MCMC: This function does MCMC sampling on fun, rejecting all samples that
%fall outside bounds, given an initial guess x0. An exit criteria may be
%provided which will terminate the function.
%
% Inputs:
%       Niter:      Number of simulations to perform. Generally should be
%                   at least 1 million.
%       stepsize:   characteristic step distance - currently is set to
%                   uniform step, but could be modified to a Gaussian.
%       fun:        function to be evaluated. Should take as inputs: first
%       x0, then bounds, then the variables in varargin in order.
%       x0:         initial guess for x (x is m x 1)
%       bounds:     a mx2 matrix specifying lower (col.1) and upper (col.2)
%                   bounds on x.
%       A:          matrix of conditions: enforced as Ax =<b
%       b:          vector of conditions: enforced as Ax =<b
%       bcut:       number of samples to discard at the beginning of
%                   sampling, known as 'burn-in'. Can use order 10% of
%                   Niter.
%       k:          Number of samples to jump in Markov chain; i.e. the
%                   algorithm will return every kth sample. This is done
%                   because consequetive samples will be correlated. I
%                   often use k = 100 if Niter = 1 million. Probably >100
%                   is not necessary, and less is probably fine in most
%                   cases.
%       exit_cond:  an exit condition that, if met, will terminate the
%                   sampler. This is in the form of a number of accepted
%                   samples. Note that you should not use this in general.
%       varargin:   input variables for fun
%
% Outputs:
%       xhats:      MCMC posterior distributions on x. Will have dimension
%                   Niter/k - bcut x m.
%       all_likes:  The likelihoods of the xhats
%       all_dpreds: Predicted data for each xhat
%       accept_rat: Acceptance ratio (# of kept samples vs. # of discarded
%                   samples)

% Created by Jeremy Maurer, Stanford, January 2016
% Last modified by Lucile Bruhat, December 6, 2016

nprint = varargin{1};
param = varargin{2};
name = varargin{3};
flag_plots = varargin{4};

models = struct;

%% specify values
N= length(x0);

%% specify default conditions
if isempty(exit_cond)
    exit_cond = inf;
end

if isempty(x0)
    x0 = (bounds(:,2) - bounds(:,1))./2;
end

if isempty(A)
    A = eye(N);
    b = inf([N, 1]);
end

if isempty(bcut)
    bcut = 10;
end


%% intialize loop
xprev = x0;
[loglikeprev, dpred,~] = feval(fun, x0, bounds, varargin{:});
lb = bounds(:,1);
ub = bounds(:,2);

xhats = nan(N, Niter/k);
%all_dpreds = nan(length(dpred), Niter/k);
[accept_rats, all_likes] = deal(nan(Niter/k, 1));

accept_count = 0; % sample accepted after Markov Chain selection
loop_ok = 0; % total number of samples that didn't fail the forward model (can be accpeted or not within the MCMC afterward) 

if fail == 1;
    % Keep track of  forward model failure
    modfails = zeros(length(x0)+1,2);
    fails = 0;
end

%% Run MCMC loop
for loop=1:Niter
    
    x = xprev + stepsize.*(rand(size(xprev))-.5);
    
    mins = x<lb;
    maxs = x>ub;
    testbd = A*x;
    test = testbd > b;

    % test if bounds are violated
    if any(mins) || any(maxs)
        accept = 0;
        
        % test secondary criteria
    elseif any(test)
        accept = 0;
    else
        [loglike, dpred,flag] = feval(fun, x, bounds, varargin{:});
        lograt = exp(loglike - loglikeprev);

        if lograt>1
            accept=1;
        else
            r=rand;
            if r<lograt
                accept=1;
            else
                accept=0;
            end
        end
        
        % Don't want to accept model when forward calculation failed
        if (flag ~= 1)
            % save xprop, flag
            if fail == 1;
                fails = fails+1;           
                modfails(:,fails) = [flag; x];
            end
            accept = 2;
        end
    end
    
    if accept==1; % update accept_count & loop_ok
        xprev = x;
        loglikeprev = loglike;
        accept_count = accept_count+1;
        loop_ok = loop_ok+1;
    elseif accept==0 % update only loop_ok
        x = xprev;
        loglike = loglikeprev;
        loop_ok = loop_ok + 1;
    elseif accept==2 % failure of the forward model, try again.
        x = xprev;
        loglike = loglikeprev;
    end
    
    %save every kth sample
    if ((mod(loop_ok,k)==0)&&(loop_ok~=0))
        xhats(:,loop_ok/k) = x;
        all_likes(loop_ok/k) = loglike;
        %all_dpreds(:,loop_ok/k) = dpred;
        accept_rats(loop_ok/k) = accept_count/loop_ok;
    end
    
    % exit_cond command section: can edit as needed
    if accept_count > exit_cond
        break
    end
    
    if (mod(loop,nprint) == 0)
        display(' ')
        display([num2str(loop/Niter*100),'% (',num2str(loop),'/',num2str(Niter),')']);        
        
        %save every nprint iterations (usually 1000)
        models.x = xhats;
        models.L = all_likes;
        %models.dpred = all_dpreds;
        models.cnt = accept_rats;
        models.xbnds = bounds;
        models.xstep = stepsize;
        if fail == 1;
        models.fails = modfails;
        end
        models.param = param;
        
        save(['models_inprogress',name,'.mat'], '-struct', 'models')
        llk = all_likes(~isnan(all_likes));
        display(['Current acceptance rate = ', num2str(100*accept_count/loop_ok,2),'%'])
        display(['Samples that failed the forward model = ', num2str((1-loop_ok/loop)*100),'%'])

        if ~isempty(llk)
            display(['Current loglike =',num2str(llk(end)),...
                '& max loglike =',num2str(max(models.L))]);
        else
            display('loglike is empty, still looking for a solution')
        end
        
    end
    
    
end

%% Remove nan is xhats
[~, col] = find(isnan(xhats));
xhats(:,unique(col))=[];
all_likes(unique(col))=[];
%all_dpreds(:,unique(col))=[];
display(' ')
display(['Samples that failed the forward model = ', num2str((1-loop_ok/loop)*100),'%'])


%% Trim initial burn-in period
xhats(:,1:bcut) = [];
all_likes(1:bcut) = [];
%all_dpreds(:,1:bcut) = [];
accept_rats(1:bcut) = [];
%accept_rat = accept_rats(end);

% save the results of the MCMC in a structure models

models.x = xhats;
models.L = all_likes;
%models.dpred = all_dpreds;
models.cnt = accept_rats;
models.xbnds = bounds;
models.xstep = stepsize;
if fail == 1;
models.fails = modfails;
end
models.param = param;

save(['models_final',name,'.mat'], '-struct', 'models')
display(' ')
end
