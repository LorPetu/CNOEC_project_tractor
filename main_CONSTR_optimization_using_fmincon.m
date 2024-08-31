%% Initialization
clear all
close all
clc
warning('off', 'all');
addpath(genpath('Models/'));
addpath(genpath('Optimization_functions/'));


error=0;
%% Model Parameters

Lt        =   3;                  % Wheelbase (m)
Li        =   2;                  % Wheelbase of implements
d         =   4.0;                  % Row width (m)

parameters=[Lt;Li;d];

%% Mode selection
% '00' - Only tractor model
% '01' - Tractor and implement model

MODE    = '01';

%% Boundaries

% Upper bound y<mx+q
constr_param.m(1)   =  0; % zero for standard case
constr_param.q(1)   = 20;

% Lower bound y<mx+q
constr_param.m(2)   =   0; % zero for standard case
constr_param.q(2)   =   0;



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

vsat        =   15/3.6;              % Input saturation
asat        =   1;                   % Cart position limits
deltasat    =   30*pi/180;
delta_psi_sat = 75*pi/180;

tol_f = [0.05,0.05,5*pi/180,0.5/3.6,0.05,0.05,5*pi/180,0.5/3.6]'; % Tolerances for the final state error

constr_param.vsat           =   vsat;
constr_param.delta_psi_sat  =   delta_psi_sat;
constr_param.tol_f          =   tol_f;



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
    'MaxFunctionEvaluations',1e6, ...
    'MaxIterations',200,...
    'StepTolerance',1e-8,...
    'HessianApproximation', 'bfgs', ...
    'PlotFcn', {@plotfun_tractor_states},... 
    'Display','iter');

%% Run solver
tic ;

constr_param.c_vel = 0; 

 U0              = [0.5*ones(8,1);
                   -0.5*ones(Np-8,1); 
                   0.2*ones(ceil(Np/2),1);
                   -0.2*ones(floor(Np/2),1);
                   Ts;]; 

[z01] = Tractor_traj(U0,z0,Nu,Ns,parameters,MODE);


[Ustar,fxstar,niter,exitflag,xsequence] = fmincon(@(U)cost_tractor_mincon(U,z0,Nu,Ns,parameters,constr_param,MODE)...
                                                    ,U0,[],[],[],[],lb,ub,...
                                                    @(U)constr_tractor_mincon(U,z0,Nu,Ns,parameters,constr_param,MODE),options);

disp(['Vincolo sul limite superiore è ', num2str(Ustar(end-1)) ]);
%% eventuale seconda iterazione
s =    Ustar(2*Np+1:end-1,1);
if exitflag.constrviolation >options.ConstraintTolerance
     constr_param.c_vel = 1; 

    U0              = [-0.5*ones(12,1);
                    0.5*ones(12,1);
                   -0.5*ones(Np-24,1); 
                   0.2*ones(10,1);
                   -0.7*ones(Np-27,1);
                   0.5*ones(17,1);
                   Ts;]; 

    
    [z02] = Tractor_traj(U0,z0,Nu,Ns,parameters,MODE);

    [Ustar,fxstar,niter,exitflag,xsequence] = fmincon(@(U)cost_tractor_mincon(U,z0,Nu,Ns,parameters,constr_param,MODE)...
                                                    ,U0,[],[],[],[],lb,ub,...
                                                    @(U)constr_tractor_mincon(U,z0,Nu,Ns,parameters,constr_param,MODE),options);

    if exitflag.constrviolation >options.ConstraintTolerance
        error=1;
    end
end

tempo_trascorso = toc;

%% calcolo stati finali
if error==0
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
    
    
    
    
    close all
    
    figure(3)
    subplot(4,1,1),plot(0:Ts:Ns*Ts,ang,'b','linewidth',1),xlabel('Time [s]'),ylabel('tractor heading angle [rad]', 'Rotation', 0),grid on
    subplot(4,1,2),plot(0:Ts:Ns*Ts,vel,'b','linewidth',1),xlabel('Time [s]'),ylabel('tractor speed [m/s]', 'Rotation', 0),grid on
    sgtitle('Heading anlge and speed of tractor and implement')
   
    figure(2)
    
    subplot(2,1,1);plot(0:Ts_p:(Np-1)*Ts_p,delta,'b','linewidth',1),xlabel('Time [s]'),ylabel('steering angle  [rad]', 'Rotation', 0),grid on
    subplot(2,1,2);plot(0:Ts_p:(Np-1)*Ts_p,acc,'b','linewidth',1),xlabel('Time [s]'),ylabel('acceleration [m/s^2]', 'Rotation', 0),grid on;
    sgtitle('Input variables')

    figure(1)
    plot(plx,ply,'Color','b','DisplayName', 'Tractor trajectory', 'linewidth',1);hold on;
    plot(asse,constr_param.m(2)*asse + constr_param.q(2),"r",'DisplayName','Lower limit' ,'linewidth',1); hold on;
    plot(asse,constr_param.m(1)*asse + constr_param.q(1),"r",'DisplayName', 'Upper limit','linewidth',1); hold on;
    daspect([1 1 1]);axis([-5 10 -5 constr_param.q(1)+5]);
    xlabel('X [m]'); ylabel('Y [m]', 'Rotation', 0);sgtitle('Trajectory'),grid on
    
    xLimits = get(gca, 'XLim');
    yLimits = get(gca, 'YLim');
    set(gca, 'XTick', xLimits(1):5:xLimits(2));
    set(gca, 'YTick', yLimits(1):5:yLimits(2))
    legend('show');
    
    
    if strcmp(MODE,'00')
        figure(1)
        plot(zf(1),zf(2),"xr",'MarkerSize', 10, 'LineWidth', 2,'DisplayName', 'Target point');
        legend('show');
    
    end
    
    if strcmp(MODE,'01')
        delta_psi=abs(zstar(3,:)-zstar(7,:));
    
        figure(4)
        plot(0:Ts:(Ns)*Ts,delta_psi,'linewidth',1,'DisplayName', 'Relative angle'); hold on
        plot(0:Ts:(Ns)*Ts,delta_psi_sat+0*(0:Ts:(Ns)*Ts),'linewidth',1,'DisplayName', 'Maximum error')
        xlabel('Time [s]'); ylabel('Relative angle [rad]', 'Rotation', 0);sgtitle('Relative angle between tractor and implement'),grid on
        legend('show');

        figure(3)
        subplot(4,1,3),plot(0:Ts:Ns*Ts,zstar(7,:),'g','linewidth',1),xlabel('Time [s]'),ylabel('implement heading angle [rad] ', 'Rotation', 0),grid on
        subplot(4,1,4),plot(0:Ts:Ns*Ts,zstar(8,:),'g','linewidth',1),xlabel('Time [s]'),ylabel('implement speed [m/s]', 'Rotation', 0),grid on
    
        figure(2)
        
        figure(1)
        plot(zstar(5,:),zstar(6,:),'Color','g','DisplayName', 'Implement', 'linewidth',1); hold on
        plot(zf(5),zf(6),"xr",'MarkerSize', 10, 'LineWidth', 2,'DisplayName', 'Target point');
        legend('show');
    
    end
    
    
    fprintf('\n    xt finale    xt target     error\n %f    %f    %f\n',plx(end,1),zf(1),abs(plx(end,1)-zf(1)));
    fprintf('   yt finale    yt target     error\n %f    %f    %f\n',ply(end,1),zf(2),abs(ply(end,1)-zf(2)));
    fprintf('   psit finale  psit target   error\n %f    %f    %f (=%f deg)\n',ang(end,1),psitf,abs(ang(end,1)-psitf),abs(ang(end,1)-psitf)*180/pi);
    fprintf('   vt finale    vt target     error\n %f    %f    %f\n\n',vel(end,1),vtf,abs(vel(end,1)-vtf));
    if MODE=='01'
        fprintf('   xi finale    xi target     error\n %f    %f    %f\n',zstar(5,end),zf(5),abs(zstar(5,end)-zf(5)));
        fprintf('   yi finale    yi target     error\n %f    %f    %f\n',zstar(6,end),zf(6),abs(zstar(6,end)-zf(6)));
        fprintf('   psii finale  psii target   error\n %f    %f    %f (=%f deg)\n',zstar(7,end),psiif,abs(zstar(7,end)-psiif),abs(zstar(7,end)-psiif)*180/pi);
        fprintf('   vi finale    vi target     error\n %f    %f    %f\n\n',zstar(8,end),vif,abs(zstar(8,end)-vif));
    end



%% errore feasibility

elseif error==1
    clc
    close all
    disp('No feasible solution has been found')
end





