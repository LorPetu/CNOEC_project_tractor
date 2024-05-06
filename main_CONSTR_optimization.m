%% Initialization
clear all
close all
clc

i=sqrt(-1);

addpath(genpath('Models/'));
addpath(genpath('Optimization_functions/'));

% Load of the precomputed 

%load('Ustar.mat');
%% Model Parameters

Lt        =   1.85;                 % Wheelbase (m)
Hti       =   pi/2;                 % Initial heading of the tractor (rad)
Htf       =   pi/2;                 % Final heading of the tractor (rad)
d         =   5.0;                  % Row width (m)
Li        =   2.5;                  % Wheelbase of implements


parameters=[Lt;Hti;Htf;d;Li];

%% initial states
xt        =      0;                 % inertial X position (m)
yt        =      0;                 % inertial Y position (m)
psit      =      pi/2;              % yaw angle (rad)
vt        =      4;                     % body x velocity (m/s) 

z0=[xt;yt;psit;vt];

%% final state
xf      =       3;
yf      =       0;
psif    =       -pi/2;
vf      =       4;

zf      =       [xf;yf;psif;vf];

%% Control problem parameters
Ts_p          =   0.25;                       % Sampling time of computation of input
% Ts_s         =   0.025;                        % Sampling time of the simulation

Tend        =   4;                         % Time horizon
% Ns          =   Tend/Ts_s;                    % Simulation steps
Np          =   Tend/Ts_p;                   % Prediction steps

time_FFD    =   [0:0.01:(Np-1)*Ts_p];
Nblock      =   Ts_p/0.01;
Nsim_FFD    =   length(time_FFD);

vsat        =   10;                     % Input saturation
asat        =   3;                      % Cart position limits
deltasat    =   30*pi/180;

% Optimization_opt

Optimization_opt.Ts_p   = Ts_p;
Optimization_opt.Tend   = Tend;

%% Initial guess

 U0          = [ zeros(Np,1);       % Accelleration (m/s^2),
                zeros(Np,1)];      
%U0=Ustar

%% Linear Constraints

C       =       [-eye(2*Np)
                eye(2*Np)];
d       =       [-deltasat*ones(Np,1);
                 -asat*ones(Np,1);
                 deltasat*ones(Np,1);
                 asat*ones(Np,1);];

% Number of equality constraints
p = 4;

% Number of inequality constraints
q = 2*Nsim_FFD;

%% Solution -  BFGS
% Initialize solver options
myoptions               =   myoptimset_const;
myoptions.Hessmethod  	=	'BFGS';
myoptions.gradmethod  	=	'CD';
myoptions.graddx        =	2^-17;
myoptions.tolgrad    	=	1e-8;
myoptions.ls_beta       =	0.2;
myoptions.ls_c          =	.01;
myoptions.ls_nitermax   =	1e2;
myoptions.nitermax      =	1e3;
myoptions.xsequence     =	'off';
% myoptions.outputfcn     =   @(U)Tractor_traj(U,z0,zf,Np,Ns,parameters,Optimization_opt);

% Run solver
[xstar,fxstar,niter,exitflag,xsequence] = myfmincon(@(U)tractor_cost_constr(U,Ts_p,Np,parameters),U0,[],[],C,d,p,q,myoptions);


