% MCMC inversion for the surface deformation caused earthquake cycle
% model rupturing down to depth D and creeping between D and H
%
% Code used in Bruhat 
%
% Lucile Bruhat
% Created: August 2018
% Last modified: August 2019

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

% Running time: XXmin for XX iterations

close all;

% Initialization
clc;
rng(1989)
addpath('mcmc/')
addpath('fracturemechanics/')
addpath('crackmodeling/')
addpath('viscoelastic/')
addpath('plot/')

dataname = 'Carrizo';

% Simulation name
name = ['_propagatingcreep_',char(dataname),'_08022019'];
display(name(2:end))

% Define structure that will transporte all the useful parameters
param = struct;

% set timer
tic;

%% Data to invert

load('data/SCEC_Carrizo_data_3.mat'); % load data
load('data/corr_hosgri.mat'); % load Hosgri fault correction

% Account for 3D correction
correct3D = 1;

if correct3D == 1
    load('data/3Deffects_3.mat'); % load corection for 3D effects
    dobs = dsc - corr + corr_hosgri;    % Velocities (mm/yr)
else
    dobs = dsc + corr_hosgri;           % Velocities (mm/yr)
end

xobs = xsc;          % Distances normal to fault (km)
error = sigsc;         % Uncertainties in the velocities (mm/yr)
Cd = diag(error.^2);               % Covariance matrix
nobs = length(xobs);

% Store in param
param.dobs = dobs;
param.Cd = Cd;
param.xobs = xobs;
param.nobs = nobs;
param.error = error;

% time since last earthquake
t = 161;
param.t = t;

param.tol = 0.05;

%% Parameters for the crack modeling

% Numerical parameters
M = 200; % number of points in z
N = 1; % number of considered Tchebyshev polynomials
mu = 3e10; % shear modulus in Pa
nu = 0.25;

% antiplane vs in-plane deformation
param.mode = 3; % antiplane
if (param.mode == 2)
    mu_s = mu/(1-nu);
elseif (param.mode == 3)
    mu_s = mu;
else
    warning('Unknown fracture mode')
end

% cn and cndot are zeros (for i>1)
cn = zeros(N-1,1);
cndot = zeros(N-1,1);

% Store in param
param.M = M;
param.N = N;
param.mu = mu;
param.nu = nu;
param.cn = cn;
param.cndot = cndot;

%% Green's functions

% Elastic (for creeping section)
Gc = zeros(nobs,M-1);

% z coordinates in m
z = linspace(2e3,50e3,M)';

% Compute Green's functions between surface rates & slip rates
for i = 1:M-2
    Gc(:,i) = atan(xobs./(z(i+1)*1e-3))/pi - atan(xobs./(z(i)*1e-3))/pi;
end
Gc(:,end) = - atan(xobs./(z(end)*1e-3))/pi;

param.Gc = Gc;
param.z = z;

%% Stress Green's function

dip = 90;
[Gs, ~] = BEM_thrust(z/sind(dip), dip, mu, nu);
param.Gs = Gs;

%% set the MCMC parameters

% Number of simulations to perform
Niter= 1e3;
% Print every nprint iteration
nprint = 100;
% Function to be evaluated
fun =  'predict_data_propagating_creep';

%  initial guess for x (x is m x 1)
x0 = [9.9683   18.3337  33.8   7.9771   70.1288   10.0041   17.3577   29.0365]';
% Characteristic step distance
stepsize = 0.3*[1   1   1  1  5   1   1   1]';
%  bounds (lower (col.1) and upper (col.2)) on x
bounds = [5 10;18 100;29 37;4 8;0 500;5 15;15 50;-50 50];
% Conditional matrix, not used here
A = [];b=[];
% Number of samples to discard at the beginning ('burn-in').
bcut = 1;
% Number of samples to jump in Markov chain (Niter/k has to be an integer)
k = 10;
% Exi; condition, not used here
exit_cond = [];
fail = 0;

%% Run the MCMC

display(['Starting MCMC procedure with Niter = ',num2str(Niter,2)])
models = mcmc_new(Niter, stepsize, fun, x0, bounds,A, b, bcut, k, exit_cond,fail,nprint,param,name,[],[]);

%% Results

% Model with highest likelihood
x_maxlike = models.x(:,models.L == max(models.L));
x_maxlike = x_maxlike(:,1);

models.x(9,:) = (models.x(6,:)-models.x(1,:))*1e3./(models.x(4,:)/(models.x(3,:)/1000)-161);
x_maxlike(9) = (x_maxlike(6)-x_maxlike(1))*1e3/(x_maxlike(4)/(x_maxlike(3)/1000)-161);

display(' ')
display('Model with highest likelihood')
display(['D = ',num2str(x_maxlike(1),3),'km'])
display(['H = ',num2str(x_maxlike(2),3),'km'])
display(['Interseismic long-term rate = ',num2str(x_maxlike(3)',3),'mm/yr'])
display(['EQ coseismic displacement = ',num2str(x_maxlike(4),3),'m'])
display(['Relaxation time = ',num2str(x_maxlike(5),3),'yrs'])
display(['Locking depth = ',num2str(x_maxlike(6),3),'km'])
display(['Depth of Hcreep = ',num2str(x_maxlike(7),3),'km'])
display(['Propagation velocity = ',num2str(x_maxlike(9),3),'m/yr'])

% Model with median values
[~, col] = find(isnan(models.x));
models.x(:,unique(col))=[];
x_med = median(models.x,2);

%% Example plots

ParametersName = {'Coseismic depth (km)','H (km)','vinf (mm/yr)','Coseismic slip (m)','Relaxation time (yrs)',...
    'Locking depth (km)','Hcreep','','vup'};

% Correlations
k = 10;i1 = 1e3;
N = length([1 2 3 4 5 6 7 9]);
figure
id = 1;
for i=[1 2 3 4 5 6 7 9]
    jd = 1;
    for j=[1 2 3 4 5 6 7 9]
        subplot(N,N,N*(id-1)+jd)
        if id==jd
            hist(models.x(i,1:k:end));
            hold all
            plot(x_maxlike(i)*[1 1],[0 1.1*max( hist(models.x(i,1:k:end)))],'r')
            title(ParametersName{i})
        elseif id < jd
            dscatter(models.x(j,1:k:end)',models.x(i,1:k:end)')
        end
        jd = jd+1;
    end
    id = id+1;
end

% Likelihood
figure
plot(models.L);title('log(likehood)')

% Plot model with highest likelihood

