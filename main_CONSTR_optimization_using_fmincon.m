%% Initialization
clear all
close all
clc
warning('off', 'all');
addpath(genpath('Models/'));
addpath(genpath('Optimization_functions/'));


%% Model Parameters

Lt        =   3;                  % Wheelbase (m)
Li        =   2;                  % Wheelbase of implements
d         =   4.0;                  % Row width (m)

parameters=[Lt;Li;d];

%% Boundaries

% Upper bound y<mx+q
constr_param.m(1)   =  0.5; % zero for standard case
constr_param.q(1)   = 8;

% Lower bound y<mx+q
constr_param.m(2)   =   0.5; % zero for standard case
constr_param.q(2)   =   0;

%% Mode selection
% '00' - Only tractor model
% '01' - Tractor and implement model

MODE    = '00';


%% initial states
% Tractor
xt      =  0;                   % inertial X position (m)
yt      =  0;                   % inertial Y position (m)
psit    =    pi/2;              % yaw angle (rad)
vt      =    4/3.6;             % body x velocity (m/s) 

% Implement 
xi      =   0;                  % implement inertial X position (m)
yi      =   0;                  % implement inertial Y position (m)
psii    =   psit;               % implement yaw angle (rad)
vi      =   vt;                 % implement body x velocity (m/s)

z0      =   [xt;yt;psit;vt];

if strcmp(MODE,'01')
 
    z0=[xi+Li*cos(psii);yi+Li*sin(psii);psit;vt;xi;yi;psii;vi];
end



%% final state
% Tractor
xtf     =   xt + d;                                         % inertial X position (m)
ytf     =   constr_param.m(2)*xtf + constr_param.q(2);      % inertial Y position (m)
psitf   =   -pi/2;                                          % yaw angle (rad)
vtf     =   4/3.6;                                          % body x velocity (m/s) 

% Implement
xif     =    xi+d;                                          % implement inertial X position (m)
yif     =    constr_param.m(2)*xif + constr_param.q(2);     % implement inertial Y position (m)
psiif   =    psitf;                                         % implement yaw angle (rad)
vif     =    vt;                                            % implement body x velocity (m/s)

zf      =    [xtf;ytf;psitf;vtf];

if strcmp(MODE,'01')

     zf=[xif+Li*cos(psiif); yif+Li*sin(psiif); psitf;vtf;xif;yif;psiif;vif];
end

constr_param.zf = zf;
%% Control problem parameters

Ns          =   75;                  % Simulation steps
Ts          =   0.25;                % initial guess for time step
Nu          =   2;                   %ogni quanti istanti di simulazione viene calcolato u

vsat        =   20/3.6;              % Input saturation
asat        =   1;                   % Cart position limits
deltasat    =   30*pi/180;
delta_psi_sat = 90*pi/180;

Optimization_opt.vsat       = vsat;
Optimization_opt.deltasat   = deltasat;
Optimization_opt.asat       = asat;
Optimization_opt.delta_psi_sat=delta_psi_sat;
Optimization_opt.Ns   = Ns;
Optimization_opt.Nu   = Nu;

Np=ceil((Ns+1)/Nu);


%% Linear Constraints

lb       =       [-deltasat*ones(Np,1);
                 -asat*ones(Np,1)];

ub        =        [deltasat*ones(Np,1);
                   asat*ones(Np,1)];
 

%% Matlab fmincon options

options = optimoptions(@fmincon,...
    'Algorithm','interior-point',...
    'FiniteDifferenceType','central',...
    'ConstraintTolerance', 1e-4,... 
    'FunctionTolerance',1e-12,...
    'EnableFeasibilityMode', true,...
    'MaxFunctionEvaluations',1e10, ...
    'MaxIterations',500,...
    'StepTolerance',1e-17,...
    'OptimalityTolerance',1e-20,...  
    'HessianApproximation', 'bfgs', ...
    'PlotFcn', {@plotfun_tractor_states},... 
    'Display','iter-detailed');

%% Run solver
tic ;

constr_param.lb_vel = 0; 

 U0              = [0.5*ones(8,1);
                   -0.5*ones(Np-8,1); 
                   0.2*ones(floor(Np/2),1);
                   -0.2*ones(ceil(Np/2),1);
                   Ts;]; 


[Ustar,fxstar,niter,exitflag,xsequence] = fmincon(@(U)cost_tractor_mincon(U,z0,parameters,Optimization_opt,constr_param,MODE)...
                                                    ,U0,[],[],[],[],lb,ub,...
                                                    @(U)constr_tractor_mincon(U,z0,parameters,Optimization_opt,constr_param,MODE),options);

disp(['Vincolo sul limite superiore è ', num2str(Ustar(end-1)) ]);
%% eventuale seconda iterazione
s =    Ustar(2*Np+1:end-1,1);
if exitflag.constrviolation >options.ConstraintTolerance
     constr_param.lb_vel = 1; 

    U0              = [-0.5*ones(12,1);
                    0.5*ones(12,1);
                   -0.5*ones(Np-24,1); 
                   0.2*ones(10,1);
                   -0.7*ones(Np-27,1);
                   0.5*ones(17,1);
                   Ts;]; 

    [Ustar,fxstar,niter,exitflag,xsequence] = fmincon(@(U)cost_tractor_mincon(U,z0,parameters,Optimization_opt,constr_param,MODE)...
                                                    ,U0,[],[],[],[],lb,ub,...
                                                    @(U)constr_tractor_mincon(U,z0,parameters,Optimization_opt,constr_param,MODE),options);

    disp(['Vincolo sul limite superiore è ', num2str(Ustar(end-1)) ]);
end

tempo_trascorso = toc;

%% calcolo stati finali

[zstar] = Tractor_traj(Ustar,z0,Nu,Ns,parameters,MODE);
 
Ts_p= Ustar(end,1)*Nu;
Ts  = Ustar(end,1);

% Visualizza il tempo trascorso
disp(['Tempo finale Tend: ', num2str(Ts*Ns), ' secondi']);
disp(['Tempo di campionamento Ts: ', num2str(Ts), ' secondi']);
disp(['Tempo per calcolo: ', num2str(tempo_trascorso), ' secondi']);
plx =   zstar(1,:)';
ply =   zstar(2,:)';
ang =   zstar(3,:)';
vel =   zstar(4,:)';

delta   =   Ustar(1:Np,1);
acc     =   Ustar(Np+1:end-1,1);

asse=linspace(-5,10,2);






figure(2)
subplot(4,1,1),plot(0:Ts:Ns*Ts,ang,'b'),xlabel('Time (s)'),ylabel('psi tractor'),grid on
subplot(4,1,2),plot(0:Ts:Ns*Ts,vel,'b'),xlabel('Time (s)'),ylabel('velocità tractor'),grid on

figure(4)
subplot(2,1,1);plot(0:Ts_p:(Np-1)*Ts_p,delta,'b'),xlabel('Time (s)'),ylabel('delta'),grid on
subplot(2,1,2);plot(0:Ts_p:(Np-1)*Ts_p,acc,'b'),xlabel('Time (s)'),ylabel('acc'),grid on;

figure(3)
plot(plx,ply,'b','DisplayName', 'Tractor');hold on;
plot(asse,constr_param.m(2)*asse + constr_param.q(2),"red",'DisplayName', 'Upper limit'); hold on;
plot(asse,constr_param.m(1)*asse + constr_param.q(1),"red",'DisplayName', 'Lower limit'); hold on;

daspect([1 1 1]);%axis([-5 10 -5 10]);
xlabel('x'); ylabel('y');title('traiettoria'),grid on
legend('show');

if strcmp(MODE,'01')
    delta_psi=abs(zstar(3,:)-zstar(7,:));

    figure(5)
    plot(0:Ts:(Ns)*Ts,delta_psi)

    figure(2)
    subplot(4,1,3),plot(0:Ts:Ns*Ts,zstar(7,:),'g'),xlabel('Time (s)'),ylabel('psi implement'),grid on
    subplot(4,1,4),plot(0:Ts:Ns*Ts,zstar(8,:),'g'),xlabel('Time (s)'),ylabel('velocità implement'),grid on

    figure(3)
    plot(zstar(5,:),zstar(6,:),'g','DisplayName', 'Implement'); hold on
    plot(zf(5),zf(6),"xr",'MarkerSize', 10, 'LineWidth', 2,'DisplayName', 'Target point');


end


% Annotation for parameters 

ann1str = sprintf('Execution time \n T_{end} = %.2f sec ',Ts*Ns); % annotation text
ann1pos = [0.018 0.71 0.22 0.16]; % annotation position in figure coordinates
ha1 = annotation('textbox',ann1pos,'string',ann1str);
ha1.HorizontalAlignment = 'left';


%Annotation for constraints
ann2str = sprintf('Constraints:\n Y < %.1f*X + %.f \n Y > %.1f*X + %.f ',constr_param.m(1),constr_param.q(1),constr_param.m(2),constr_param.q(2)); % annotation text
ann2pos = [0.02 0.2 0.1 0.1]; % annotation position in figure coordinates
ha2 = annotation('textbox',ann2pos,'string',ann2str);
ha2.HorizontalAlignment = 'left';
ha2.EdgeColor = 'red';

% fprintf('\n    xt finale    xt target     error\n %f    %f    %f\n',plx(end,1),zf(1),abs(plx(end,1)-zf(1)));
% fprintf('   yt finale    yt target     error\n %f    %f    %f\n',ply(end,1),zf(2),abs(ply(end,1)-zf(2)));
% fprintf('   psit finale  psit target   error\n %f    %f    %f (=%f deg)\n',ang(end,1),psitf,abs(ang(end,1)-psitf),abs(ang(end,1)-psitf)*180/pi);
% fprintf('   vt finale    vt target     error\n %f    %f    %f\n\n',vel(end,1),vtf,abs(vel(end,1)-vtf));
% fprintf('   xi finale    xi target     error\n %f    %f    %f\n',zstar(5,end),zf(5),abs(zstar(5,end)-zf(5)));
% fprintf('   yi finale    yi target     error\n %f    %f    %f\n',zstar(6,end),zf(6),abs(zstar(6,end)-zf(6)));
% fprintf('   psii finale  psii target   error\n %f    %f    %f (=%f deg)\n',zstar(7,end),psiif,abs(zstar(7,end)-psiif),abs(zstar(7,end)-psiif)*180/pi);
% fprintf('   vi finale    vi target     error\n %f    %f    %f\n\n',zstar(8,end),vif,abs(zstar(8,end)-vif));

