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
Ts_p          =   0.25;                       % Sampling time of comutation of input
Ts_s         =   0.025;                        % Sampling time of the simulation
Q           =   diag([1;1;1;1]);       % Tracking error weight
Qf          =   diag([1;1;1;1]*1e5);       % Terminal weight
Qdot        =   diag([0;0;1;1]*10);
R           =   diag([1;1]*10);

Tend        =   4;                         % Time horizon
Ns          =   Tend/Ts_s;                    % Simulation steps
Np          =   Tend/Ts_p;                   % Prediction steps

vsat        =   10;                     % Input saturation
asat        =   3;                      % Cart position limits
deltasat    =   30*pi/180;

alpha       =   1e1;               % Barrier function scaling
beta        =   1e1;                 % Barrier function coefficient

% Optimization_opt

Optimization_opt.Q          = Q;
Optimization_opt.Qf         = Qf;
Optimization_opt.Qdot       = Qdot;
Optimization_opt.R          = R;
Optimization_opt.vsat       = vsat;
Optimization_opt.deltasat   = deltasat;
Optimization_opt.asat       = asat;
Optimization_opt.alpha      = alpha;
Optimization_opt.beta       = beta;
Optimization_opt.Ts_s   = Ts_s;
Optimization_opt.Ts_p   = Ts_p;
Optimization_opt.Tend   = Tend;

%% Initial guess

 U0          = [ zeros(Np,1);       % Accelleration (m/s^2),
                zeros(Np,1)];      
%U0=Ustar

%% Solution -  BFGS
% Initialize solver options
myoptions               =   myoptimset;
myoptions.Hessmethod  	=	'BFGS';
myoptions.gradmethod  	=	'CD';
myoptions.graddx        =	2^-17;
myoptions.tolgrad    	=	1e-8;
myoptions.ls_beta       =	0.2;
myoptions.ls_c          =	.01;
myoptions.ls_nitermax   =	1e2;
myoptions.nitermax      =	1e3;
myoptions.xsequence     =	'off';
myoptions.outputfcn     =   @(U)Tractor_traj(U,z0,zf,Np,Ns,parameters,Optimization_opt);

% Run solver
[xstar,fxstar,niter,exitflag,xsequence] = myfmincon(@(x)tractor_cost_constr(x,Ts,Np,th),x0,[],[],C,d,0,q,myoptions);


