%% Initialization
clear all
close all
clc
warning('off', 'all');
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

Ns          =   50;                    % Simulation steps
Ts         =   0.25;             % initial guess for time step
Nu=2;  %ogni quanti istanti di simulazione viene calcolato u

vsat        =   20/3.6;                     % Input saturation
asat        =   1;                      % Cart position limits
deltasat    =   30*pi/180;

Optimization_opt.vsat       = vsat;
Optimization_opt.deltasat   = deltasat;
Optimization_opt.asat       = asat;

Optimization_opt.Ns   = Ns;
Optimization_opt.Nu   = Nu;

Np=ceil((Ns+1)/Nu);
s_number=4+1;

%% Linear Constraints

lb       =       [-deltasat*ones(Np,1);
                 -asat*ones(Np,1);
                 zeros(s_number+1,1)];

ub        =        [deltasat*ones(Np,1);
                   asat*ones(Np,1)];

%% Parametrized Constraint
% boundaries are expressed as y=mx+q

% Upper bound y<mx+q
constr_param.m(1)   =  0.5; % zero for standard case
constr_param.q(1)   = 6;

% Lower bound y<mx+q
constr_param.m(2)   =   0.5; % zero for standard case
constr_param.q(2)   =   0; 


zf(2) = constr_param.m(2)*zf(1) + constr_param.q(2); 
constr_param.zf = zf;

% velocity
 

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
    'ConstraintTolerance', 1e-5,... 
    'FunctionTolerance',1e-6,...
    'EnableFeasibilityMode', true,...
    'MaxFunctionEvaluations',10e5, ...
    'MaxIterations',500,...
    'StepTolerance',1e-5,...
    'OptimalityTolerance',1e-12,...    
    'PlotFcn', {@plotfun_tractor_traj,@optimplotfval},... %,@optimplotfval
    'Display','iter-detailed');

%% Run solver
tic ;

constr_param.lb_vel = 0; 

 U0              = [0.5*ones(5,1);
                   -0.5*ones(Np-5,1); 
                   0.1*ones(floor(Np/2),1);
                   -0.1*ones(ceil(Np/2),1);
                   zeros(s_number,1)
                   Ts;]; 
 


[Ustar,fxstar,niter,exitflag,xsequence] = fmincon(@(U)cost_tractor_mincon(U,z0,parameters,Optimization_opt,constr_param)...
                                                    ,U0,[],[],[],[],lb,ub,...
                                                    @(U)constr_tractor_mincon(U,z0,parameters,Optimization_opt,constr_param),options);

disp(['Vincolo sul limite superiore è ', num2str(Ustar(end-1)) ]);
s =    Ustar(2*Np+1:end-1,1);
if Ustar(end-1)>1e-2 || s(1)+s(2)+s(3)+s(4)>0.5
    constr_param.lb_vel = 1; 

    U0              = [-0.4*ones(2,1);
                   -0.4*ones(Np-2,1); 
                   0.3*ones(9,1);
                   -0.8*ones(Np-18,1);
                   0.3*ones(9,1);
                   zeros(s_number,1)
                   Ts;]; 

    
    [Ustar,fxstar,niter,exitflag,xsequence] = fmincon(@(U)cost_tractor_mincon(U,z0,parameters,Optimization_opt,constr_param)...
                                                    ,U0,[],[],[],[],lb,ub,...
                                                    @(U)constr_tractor_mincon(U,z0,parameters,Optimization_opt,constr_param),options);

    disp(['Vincolo sul limite superiore è ', num2str(Ustar(end-1)) ]);
end

tempo_trascorso = toc;

% Visualizza il tempo trascorso
disp(['Tempo per calcolo: ', num2str(tempo_trascorso), ' secondi']);

%% calcolo stati finali

[zstar] = Tractor_traj(Ustar,z0,zf,Nu,Ns,parameters,Optimization_opt);
 
Ts_p= Ustar(end,1)*Nu;
Ts  = Ustar(end,1);

% Visualizza il tempo trascorso
disp(['Tempo finale Tend: ', num2str(Ts*Ns), ' secondi']);
disp(['Tempo di campionamento Ts: ', num2str(Ts), ' secondi']);

plx =   zstar(1,:)';
ply =   zstar(2,:)';
ang =   zstar(3,:)';
vel =   zstar(4,:)';

delta   =   Ustar(1:Np,1);
acc     =   Ustar(Np+1:end-s_number-1,1);


figure(2)
subplot(4,1,1),plot(0:Ts:Ns*Ts,ang),xlabel('Time (s)'),ylabel('psi'),grid on
subplot(4,1,2),plot(0:Ts:Ns*Ts,vel),xlabel('Time (s)'),ylabel('velocità'),grid on
subplot(4,1,3);plot(0:Ts_p:(Np-1)*Ts_p,delta),xlabel('Time (s)'),ylabel('delta'),grid on
subplot(4,1,4);plot(0:Ts_p:(Np-1)*Ts_p,acc),xlabel('Time (s)'),ylabel('acc'),grid on;


figure(3)
plot(plx,ply,'o');hold on;
plot(plx,constr_param.m(2)*plx + constr_param.q(2),"red"); hold on;
plot(plx,constr_param.m(1)*plx + constr_param.q(1),"red"); hold on;
plot(zf(1),zf(2),"xr",'MarkerSize', 10, 'LineWidth', 2);daspect([1 1 1]);xlabel('x'); ylabel('y');title('traiettoria'),grid on
%

% Annotation for parameters 
ann1str = sprintf('Tempo impiegato \n T_{end} = %.2f sec ',Ts*Ns); % annotation text
ann1pos = [0.018 0.71 0.22 0.16]; % annotation position in figure coordinates
ha1 = annotation('textbox',ann1pos,'string',ann1str);
ha1.HorizontalAlignment = 'left';


%Annotation for constraints
ann2str = sprintf('Constraints:\n Y < %.1f*X + %.f \n Y > %.1f*X + %.f ',constr_param.m(1),constr_param.q(1),constr_param.m(2),constr_param.q(2)); % annotation text
ann2pos = [0.02 0.2 0.1 0.1]; % annotation position in figure coordinates
ha2 = annotation('textbox',ann2pos,'string',ann2str);
ha2.HorizontalAlignment = 'left';
ha2.EdgeColor = 'red';

 fprintf('   xfinale    xtarget     error\n %f    %f    %f\n',plx(Ns,1),zf(1),abs(plx(end,1)-zf(1)));
fprintf('   yfinale    ytarget     error\n %f    %f    %f\n',ply(Ns,1),zf(2),abs(ply(end,1)-zf(2)));
fprintf('   psifinale  psitarget   error\n %f    %f    %f\n',ang(Ns,1),psif,abs(ang(end,1)-psif));
fprintf('   vfinale    vtarget     error\n %f    %f    %f\n\n',vel(Ns,1),vf,abs(vel(end,1)-vf));


