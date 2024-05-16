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
Ts_p          =   0.5;                       % Sampling time of comutation of input
Ts_s         =   0.025;                        % Sampling time of the simulation
Q           =   diag([1;1;1;1]*1e2);       % Tracking error weight
Qf          =   diag([1;1;1;1]*1e20);       % Terminal weight
Qdot        =   diag([0;0;1;1]*0);
R           =   diag([1;1]*0);

Tend        =   5;                         % Time horizon
Ns          =   Tend/Ts_s;                    % Simulation steps
Np          =   Tend/Ts_p;                   % Prediction steps

vsat        =   10;                     % Input saturation
asat        =   3;                      % Cart position limits
deltasat    =   30*pi/180;

alpha       =   1e3;               % Barrier function scaling
beta        =   1e2;                 % Barrier function coefficient

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


%% Optimization parameters

myoptions               =   myoptimset;
myoptions.Hessmethod  	=	'GN';
myoptions.GN_funF       = @(U)tractor_cost_GN_grad_mod(U,z0,zf,Np,Ns,parameters,Optimization_opt);
myoptions.gradmethod  	=	'CD';
myoptions.graddx        =	2^-26;
myoptions.tolgrad    	=	1e-10;
myoptions.tolfun        =	1e-12;
myoptions.ls_tkmax      =   1;
myoptions.ls_beta       =	0.5;   %0.5
myoptions.ls_c          =	.1;     %0.1
myoptions.ls_nitermax   =	100;     %20
myoptions.nitermax      =	100;
myoptions.xsequence     =	'on';
myoptions.outputfcn     =   @(U)Tractor_traj(U,z0,zf,Np,Ns,parameters,Optimization_opt);


%non so perchè ma con GN Ustar sembra uscire dai limiti. va controllato se
%è solo nel grafico o è un problema proprio del calcolo


tic;
if  strcmp(myoptions.Hessmethod, 'BFGS')
 [Ustar,fxstar,k,exitflag,xsequence] = myfminunc(@(U)Tractor_cost_alternative(U,z0,zf,Np,Ns,parameters,Optimization_opt),U0,myoptions);
elseif  strcmp(myoptions.Hessmethod, 'GN')
 [Ustar,fxstar,k,exitflag,xsequence] = myfminunc(@(U)tractor_cost_star_mod(U,z0,zf,Np,Ns,parameters,Optimization_opt),U0,myoptions);
end
tempo_trascorso = toc;
%% calcolo stati finali
[zstar] = stati_finali(Ustar,z0,zf,Np,Ns,parameters,Optimization_opt);
%% 

N=length(zstar);


plx=zeros(N,1);
ply=zeros(N,1);
ang=zeros(N,1);
vel=zeros(N,1);

for j=1:1:N
   plx(j,1)=zstar(1,j);
   ply(j,1)=zstar(2,j);
   ang(j,1)=zstar(3,j);
   vel(j,1)=zstar(4,j);
end

Nu=length(Ustar)/2;

for j=1:1:Nu
   del(j,1)=Ustar(j,1);
   acc(j,1)=Ustar(Nu+j,1);
end

figure(2)
subplot(4,1,1),plot(0:Ts_s:(N-1)*Ts_s,ang),xlabel('Time (s)'),ylabel('psi'),grid on
subplot(4,1,2),plot(0:Ts_s:(N-1)*Ts_s,vel),xlabel('Time (s)'),ylabel('velocità'),grid on
subplot(4,1,3);plot(0:Ts_p:(Nu-1)*Ts_p,del),xlabel('Time (s)'),ylabel('delta'),grid on
subplot(4,1,4);plot(0:Ts_p:(Nu-1)*Ts_p,acc),xlabel('Time (s)'),ylabel('acc'),grid on;

figure(3)
plot(plx,ply,'o');hold on;
plot(xf,yf,"xr",'MarkerSize', 10, 'LineWidth', 2);daspect([1 1 1]);xlabel('x'); ylabel('y');title('traiettoria'),grid on

fprintf('   xfinale    xtarget     error\n %f    %f    %f\n',plx(N,1),xf,abs(plx(end,1)-xf));
fprintf('   yfinale    ytarget     error\n %f    %f    %f\n',ply(N,1),yf,abs(ply(end,1)-yf));
fprintf('   psifinale  psitarget   error\n %f    %f    %f\n',ang(N,1),psif,abs(ang(end,1)-psif));
fprintf('   vfinale    vtarget     error\n %f    %f    %f\n\n',vel(N,1),vf,abs(vel(end,1)-vf));



% Visualizza il tempo trascorso
disp(['Tempo trascorso: ', num2str(tempo_trascorso), ' secondi']);
