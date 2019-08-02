function [loglike, dpred,  flag_fwmodel, sfig] = predict_data_propagating_creep(m,varargin)
%
% Lucile Bruhat
% Created: November 2018
% Last modified: November 2018

% Using my code for earthquake cycle deformation

% Here we invert for:
% D = m(1), full rupture depth in km
% H = m(2), downdip limit of the transition region in km, elastic thickness
% v_inf = m(3), long-term rate in mm/yr
% C = m(4), Coseismic displacement in m,  = recurrence time T x v_inf
% tR = m(5), Relaxation time in years
% d = m(6), locking depth in km
% Hcreep = m(7), depth of uniform creep in km
% alpha = m(8), free parameter to accounf for rigid block motion, in mm/yr

%% Initialization

param = varargin{3};
flag_plots = varargin{5};
sfig = varargin{6};

% Useful variables
dobs = param.dobs;
Cd = param.Cd;
alpha = m(8); % block motion

if isempty(flag_plots)
    % Compute the average slip rate
    [vtot, flag_fwmodel, sfig] = predict_avgsliprate_propagating_creep(m,param,0,sfig);
    
    if (flag_fwmodel == 1)
        
        % Compute the predicted data
        dpred = vtot + alpha;
        % Compute log likelihood
        loglike = -0.5*((dobs-dpred)'*(Cd\(dobs-dpred)));
        
    else
        %set very low probability for unconverged model to prevent it from being accepted
        loglike = -1e10;
        dpred = nan(size(dobs));
        
    end
    
else
    % Compute the average slip rate
    [vtot, flag_fwmodel, sfig] = predict_avgsliprate_propagating_creep(m,param,1,sfig);
    
    % Compute the predicted data
    dpred = vtot + alpha;
    % Compute log likelihood
    loglike = -0.5*((dobs-dpred)'*(Cd\(dobs-dpred)));
end

end

function [vtot, flag_fwmodel, sfig] = predict_avgsliprate_propagating_creep(m,param,flag_plots, sfig)
% This function computes the surface deformation caused earthquake cycle
% model rupturing up to D and fulling creeping between D and H

% Created: November 2018
% Last modified: November 2018

% Using my code for earthquake cycle deformation

% Here we invert for:
% D = m(1), full rupture depth in km
% H = m(2), downdip limit of the transition region in km, elastic thickness
% v_inf = m(3), long-term rate in mm/yr
% C = m(4), Coseismic displacement in m,  = recurrence time T x v_inf
% tR = m(5), Relaxation time in years
% d = m(6), locking depth in km
% Hcreep = m(7), depth of uniform creep in km
% alpha = m(8), free parameter to accounf for rigid block motion, in mm/yr


%% Initialization

vtot = 0;

% Load parameters to test
D = m(1); % Full rupture depth in km
H = m(2); % Elastic thickness in km
dlock = m(6); % Current locking depth in km
Hcreep = m(7); % thickness at whcich slip rate before constant, in km

% observation points (distance in km)
x = param.xobs;

% Plate motion velocity in mm/yr
v_inf = m(3);
param.v_inf = v_inf/(1000*pi*1e7); % save in m/s

% time in yr
t = param.t;
ts = t*pi*1e7;% convert in s for crack modeling
param.t = ts;

% Relaxation parameters in yr
tR = m(5);
tRs = tR*1e7*pi; % convert in m/s
param.tR = tRs;

% Coseismic slip in m
C = m(4);

% Recurence time in yr
T = C/(v_inf/1000);
% reject if recurrence smaller than time of last earthquake
if (T<161); flag_fwmodel = 0;return;end

Ts = T*1e7*pi; % convert to seconds
param.T = Ts;

% Points on the fault
z = param.z;

% updip velocity in m/yr
vup = (dlock-D)*1e3/(T-t);

% Flags
if (D>H); flag_fwmodel = 0;return;end
if (D>dlock); flag_fwmodel = 0;return;end
if (dlock>H); flag_fwmodel = 0;return;end
if (dlock>Hcreep); flag_fwmodel = 0;return;end
if (Hcreep>H); flag_fwmodel = 0;return;end
if ((vup*t)>(H-dlock)*1e3); flag_fwmodel = 0;return;end

%% Earthquake cycle model 

% Creep model at t = T
param.t = Ts;

% Grid points for the crack

a = (Hcreep-D)*1e3; % length of the crack in m
zmax = Hcreep*1000;

if (a<=0); flag_fwmodel = 0;return;end
param.a = a;

zin = z((z>=zmax-a)&(z<=zmax)); % z inside the crack
if isempty(zin); flag_fwmodel = 0;return;end

xi = linspace(-1,1,length(zin))';

% Crack modeling

% Crack ingredients (take into account the viscoelastic relaxation)
[Db,Ddot,g,gp,fn,Fn,hn,tn] = ingredient_crack_fast_v3(xi,param);

% Store in param
param.D = Db;
param.Ddot = Ddot;
param.g = g;
param.gp = gp;
param.fn = fn;
param.Fn = Fn;
param.hn = hn;
param.tn = tn;

% cn and cndot are zeros here (first order model)
cn = 0;
cndot = zeros(size(cn));

% crack solution
[~,~,slip,~] = make_crack(xi,cn,cndot,0,param);

% Tapered coseismic distribution
zcod = zin;
slipcod = C-slip;

% Classic Savage model for viscoelastic modeling
Nn = 30;
K = 30;

zobs = z(z<H*1000);
zslip = [linspace(0,min(zcod/1000)-.5,30)';zcod/1000;linspace(max(zcod/1000)+.5,H,10)'];
coseismic = [ones(30,1)*C;slipcod;zeros(10,1)];
mu = param.mu;

% surface rates in mm/yr
vEQcycle_up = ve_cycle_cstsli(t,tR,T,Nn,K,v_inf,H,D,x,'right_lat');
vEQcycle_down = ve_cycle_varsli(t,tR,T,Nn,K,zcod/1000,slipcod,H,x,'right_lat')*1000;

vEQcycle = vEQcycle_up + vEQcycle_down;
stressrateEQcycle = se_cycle_varsli(t,tR,T,mu,Nn,K,zslip*1000,coseismic,H*1000,zobs,'right_lat')/1000;% in kPa/yr
stressrateEQcyclebottom = se_cycle_varsli_bottom(t,tR,T,mu,Nn,K,zslip*1000,coseismic,H*1000,linspace(H*1000,100*1000,50),'right_lat')/1000;% in kPa/yr

%% Deep Creep model

% Grid points for the crack

% length of the crack in m
a = (Hcreep-dlock)*1e3;
zmax = Hcreep*1000;

if (a<=0); flag_fwmodel = 0;return;end
param.a = a;

zin = z((z>=zmax-a)&(z<=zmax)); % z inside the crack
if isempty(zin); flag_fwmodel = 0;return;end

zlock = z((z<zmax-a)); % z locked region
zcreep = z((z>zmax)&(z<H*1000)); % z creeping region
zvisco = z((z>H*1000));
xi = linspace(-1,1,length(zin))';

% Crack modeling

% Crack ingredients (take into account the viscoelastic relaxation)
[Db,Ddot,g,gp,fn,Fn,hn,tn] = ingredient_crack_fast_v3(xi,param);

% Store in param
param.D = Db;
param.Ddot = Ddot;
param.g = g;
param.gp = gp;
param.fn = fn;
param.Fn = Fn;
param.hn = hn;
param.tn = tn;

% convert adot from m/yr -> m/s
adot = vup/(pi*1e7);

% Compute c2 (no stress intensity at bottom end of the crack)
cn = 0;
cndot = zeros(size(cn));

% crack solution
[~,~,~,slip_rate] = make_crack(xi,cn,cndot,adot,param);

% condition on the slip rate positiveness
if ~isempty(find(slip_rate<0,1));flag_fwmodel = 0;return;end

% condition on the slip rate maximum value (avoid rough solutions)
if ~isempty(find(slip_rate>1.1*(v_inf/(1000*pi*1e7)),1));flag_fwmodel = 0;return;end

% Condition on monotony of slip rate distribution
if max(slip_rate)>1.1*slip_rate(end);flag_fwmodel = 0;return;end

% Get average slip rate in mm/yr inside and outside the crack
avgsliprate = [zeros(size(zlock));slip_rate*1e3*1e7*pi;slip_rate(end)*1e3*1e7*pi*ones(size(zcreep));zeros(size(zvisco(2:end)))];

avgsliprate_el = [zeros(size(zlock));slip_rate*1e3*1e7*pi;slip_rate(end)*1e3*1e7*pi.*[ones(size(zcreep));ones(size(zvisco(2:end)))]];

Gs = param.Gs;
tpred = Gs(2:end,2:end)*avgsliprate_el/1000/1000; % in kPa/yr
tpredH = tpred(z<H*1000);

avgslipratein = slip_rate*1e3*1e7*pi;
%% Elastic surface rates in mm/yr

Gc = param.Gc;
velcreep = Gc*avgsliprate;

%% Viscoelastic surface rates in mm/yr
t = param.t/(pi*1e7);
xobs = x';
d = dlock-T*vup/1000;

[vvecreep,~] = ve_flow_varcreep_tvec_cycle(t,tR,Nn,zin/1000,avgslipratein,v_inf,H,D,xobs,'right_lat');
%stemp = se_flow_varcreep_tvec_cycle(t,tR,T,mu,Nn,K,zin/1000,avgslipratein,v_inf,H,D,zobs,'right_lat')/1000;

%% Sum them all
vtot = vEQcycle + velcreep + vvecreep ;
stot = stressrateEQcycle+tpredH ;
stotbottom = stressrateEQcyclebottom;

% condition on location of maximum stressing rate
if (((z(stot == max(stot))/1000)>11)&&((z(tpred == max(tpred))/1000)<7));flag_fwmodel = 0;return;end

%% Figures

if (flag_plots == 1)
    

    %% Fit to data
    dobs = param.dobs;
    xobs = param.xobs;
    error = param.error;
    alpha = m(8);
    
    % Compute the predicted data
    dpred = vtot + alpha;
    
    % save into sfig
    sfig.xobs = xobs';
    sfig.dobs = dobs';
    sfig.error = error;
    sfig.dpred = dpred';

    sfig.H = H*[1 1];
    sfig.D = D*[1 1];
    sfig.dlock = dlock*[1 1];
    
    sfig.zcod = [0 D zcod'*1e-3 H];
    sfig.slipcod = [[1 1]*C slipcod' 0];
    sfig.zsr = [0 zin'*1e-3,Hcreep-0.5,H];
    sfig.sr = [0,avgslipratein',[1 1]*avgslipratein(end)];
    sfig.zstr = [zobs'/1000,linspace(H,100,49)];
    sfig.str = [stot',stotbottom'];

    
end

flag_fwmodel = 1;

end
