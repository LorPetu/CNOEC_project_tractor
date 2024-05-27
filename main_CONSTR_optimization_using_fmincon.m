%% Initialization
clear all
close all
clc

addpath(genpath('Models/'));
addpath(genpath('Optimization_functions/'));

% Load of the precomputed 


%% Model Parameters

Lt        =   1.85;                  % Wheelbase (m)
Hti       =   pi/2;                 % Initial heading of the tractor (rad)
Htf       =   pi/2;                 % Final heading of the tractor (rad)
d         =   4.0;                  % Row width (m)
Li        =   2.5;                  % Wheelbase of implements


parameters=[Lt;Hti;Htf;d;Li];

%% initial states
xt        =      0;                 % inertial X position (m)
yt        =      0;                 % inertial Y position (m)
psit      =      pi/2;              % yaw angle (rad)
vt        =      4/3.6;                     % body x velocity (m/s) 

z0=[xt;yt;psit;vt];

%% final state
xf      =       xt+d; 
yf      =       0;
psif    =       -pi/2;
vf      =       4/3.6;

zf      =       [xf;yf;psif;vf];

%% Control problem parameters
Ts_p          =   0.3;                       % Sampling time of comutation of input
Ts_s         =   0.2; 

Tend        =   14;% NB sbagliavamo e mettevamo una                     % Time horizon

Ns          =   Tend/Ts_s;                    % Simulation steps
Np          =   ceil(Tend/Ts_p);                  % Prediction steps
vsat        =   20/3.6;                     % Input saturation
asat        =   1;                      % Cart position limits
deltasat    =   30*pi/180;

Optimization_opt.vsat       = vsat;
Optimization_opt.deltasat   = deltasat;
Optimization_opt.asat       = asat;
Optimization_opt.Ts_s   = Ts_s;
Optimization_opt.Ts_p   = Ts_p;
Optimization_opt.Tend   = Tend;
Optimization_opt.Np   = Np;
Optimization_opt.Ns   = Ns;
%% Initial guess
% U0 = load('best_initial_cond.mat').Ustar;
% U0          = best_initial_cond;
s_number=4+1;
% 
U0              = [-0.5*ones(Np,1);     %delta
                   zeros(ceil(Np/2),1);
                   -0.2*ones(floor(Np/2),1);
                   zeros(s_number,1)];   
% U0          = [-1;zeros(ceil(Np/2)-1,1);0.5;zeros(floor(Np/2)-1,1);      
%              0;zeros(ceil(Np/2)-1,1);-0.5;zeros(floor(Np/2)-1,1);
%              zeros(s_number,1)];
%              zeros(4,1)];      
% % %U0=Ustar

%% Linear Constraints

lb       =       [-deltasat*ones(Np,1);
                 -asat*ones(Np,1);
                 zeros(s_number,1)];

ub        =        [deltasat*ones(Np,1);
                   asat*ones(Np,1)];

%% Parametrized Constraint
% boundaries are expressed as y=mx+q

% Upper bound y<mx+q
constr_param.m(1)   =  0; % zero for standard case
constr_param.q(1)   = 10;

% Lower bound y<mx+q
constr_param.m(2)   =   0; % zero for standard case
constr_param.q(2)   =   0; 


zf(2) = constr_param.m(2)*zf(1) + constr_param.q(2); 
constr_param.zf = zf;

%% Solution -  BFGS
% Initialize solver options
myoptions               =   myoptimset_const;
myoptions.Hessmethod  	=	'BFGS';
myoptions.gradmethod  	=	'CD';
myoptions.graddx        =	2^-17;
myoptions.GN_funF       = @(U)tractor_cost_GN_grad_constr(U,z0,zf,parameters,Optimization_opt);
myoptions.tolgrad    	=	1e-12;
myoptions.tolfun    	=	1e-8;
myoptions.tolconstr     =   4e-2;
myoptions.ls_tkmax      =   1.3; %before 1
myoptions.ls_beta       =	0.8;
myoptions.ls_c          =	0.08;
myoptions.ls_nitermax   =	100;
myoptions.nitermax      =	200;
myoptions.xsequence     =	'on';
myoptions.outputfcn     =    @(U)Tractor_traj(U,z0,zf,Np,Ns,parameters,Optimization_opt);

%% Matlab fmincon options

   % Default properties:
   %                  Algorithm: 'interior-point'
   %         BarrierParamUpdate: 'monotone'
   %        ConstraintTolerance: 1.0000e-06
   %                    Display: 'final'
   %      EnableFeasibilityMode: 0
   %   FiniteDifferenceStepSize: 'sqrt(eps)'
   %       FiniteDifferenceType: 'forward'
   %       HessianApproximation: 'bfgs'
   %                 HessianFcn: []6
   %         HessianMultiplyFcn: []
   %                HonorBounds: 1
   %     MaxFunctionEvaluations: 3000
   %              MaxIterations: 1000
   %             ObjectiveLimit: -1.0000e+20
   %        OptimalityTolerance: 1.0000e-06
   %                  OutputFcn: []
   %                    PlotFcn: []
   %               ScaleProblem: 0
   %  SpecifyConstraintGradient: 0
   %   SpecifyObjectiveGradient: 0
   %              StepTolerance: 1.0000e-10
   %        SubproblemAlgorithm: 'factorization'
   %                   TypicalX: 'ones(numberOfVariables,1)'
   %                UseParallel: 0

options = optimoptions(@fmincon,...
    'Algorithm','interior-point',...
    'FiniteDifferenceType','central',...
    'ConstraintTolerance', 1e-6,... 
    'FunctionTolerance',1e-6,...
    'MaxFunctionEvaluations',10e5, ...
    'MaxIterations',500,...
    'StepTolerance',1e-10,...
    'OptimalityTolerance',1e-6,...    
    'PlotFcn', {@plotfun_tractor_traj,@optimplotfval},... %,@optimplotfval
    'Display','iter-detailed');

%% Run solver

tic ;
[Ustar,fxstar,niter,exitflag,xsequence] = fmincon(@(U)cost_tractor_mincon(U,z0,parameters,Optimization_opt,constr_param)...
                                                    ,U0,[],[],[],[],lb,ub,...
                                                    @(U)constr_tractor_mincon(U,z0,parameters,Optimization_opt,constr_param),options);
% [Ustar,fval,exitflag,output] = nested_fmincon(z0,parameters,Optimization_opt,constr_param,U0,C,d,[],[],[],[],options);

tempo_trascorso = toc;


%% calcolo stati finali
[zstar] = stati_finali(Ustar,z0,zf,Np,Ns,parameters,Optimization_opt);

N=length(zstar);

plx=zeros(N,1);
ply=zeros(N,1);
ang=zeros(N,1);
vel=zeros(N,1);
distN=zeros(N,1);

for j=1:1:N
   plx(j,1)=zstar(1,j);
   ply(j,1)=zstar(2,j);
   ang(j,1)=zstar(3,j);
   vel(j,1)=zstar(4,j);
end

for i=2:1:N
       distN =sqrt((plx(i,1)-plx(i-1,1))^2 + (ply(i,1)-ply(i-1,1))^2);
end

sum(distN);

Nu=Np;

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
plot(plx,ply,'o');hold on;axis([-5 10 -2 12])
plot(plx,constr_param.m(2)*plx + constr_param.q(2),"red"); hold on;
plot(plx,constr_param.m(1)*plx + constr_param.q(1),"red"); hold on;
plot(zf(1),zf(2),"xr",'MarkerSize', 10, 'LineWidth', 2);daspect([1 1 1]);xlabel('x'); ylabel('y');title('traiettoria'),grid on

% % Annotation for parameters 
% ann1str = sprintf('Opt param:\nLINE SEARCH\n tkmax = %.1f \n beta = %.1f \n c = %.2f ',myoptions.ls_tkmax,myoptions.ls_beta, myoptions.ls_c); % annotation text
% ann1pos = [0.018 0.71 0.19 0.22]; % annotation position in figure coordinates
% ha1 = annotation('textbox',ann1pos,'string',ann1str);
% ha1.HorizontalAlignment = 'left';
% 
% %Annotation for constraints
% ann2str = sprintf('Constraints:\n Y < %.1f*X + %.f \n Y > %.1f*X + %.f ',constr_param.m(1),constr_param.q(1),constr_param.m(2),constr_param.q(2)); % annotation text
% ann2pos = [0.02 0.2 0.1 0.1]; % annotation position in figure coordinates
% ha2 = annotation('textbox',ann2pos,'string',ann2str);
% ha2.HorizontalAlignment = 'left';
% ha2.EdgeColor = 'red';
% 
% %% Save figures
% m_string = replace(sprintf('m%.1f__q%.2f',constr_param.m(1),constr_param.q(1)),'.','');
% figname = sprintf('%.f___%s__q%.f ',exitflag,m_string);
% saveas(figure(3),[pwd '/Constr_sat01/' figname]);
% saveas(figure(2),[pwd '/Constr_sat01/' figname '_states']);
% 
% fprintf('   xfinale    xtarget     error\n %f    %f    %f\n',plx(N,1),zf(1),abs(plx(end,1)-zf(1)));
% fprintf('   yfinale    ytarget     error\n %f    %f    %f\n',ply(N,1),zf(2),abs(ply(end,1)-zf(2)));
% fprintf('   psifinale  psitarget   error\n %f    %f    %f\n',ang(N,1),psif,abs(ang(end,1)-psif));
% fprintf('   vfinale    vtarget     error\n %f    %f    %f\n\n',vel(N,1),vf,abs(vel(end,1)-vf));



% Visualizza il tempo trascorso
disp(['Tempo trascorso: ', num2str(tempo_trascorso), ' secondi']);

