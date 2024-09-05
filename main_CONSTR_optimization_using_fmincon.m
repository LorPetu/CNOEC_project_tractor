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
d         =   4.0;                % Row width (m)

parameters=[Lt;Li;d];

%% Mode selection
% '00' - Only tractor model
% '01' - Tractor and implement model

MODE    = '01';

%% Boundaries
% Linear approximation has been used to represent the operational area 
% delimited by the headland limits. Here are specified the parameters 

% Upper bound y<mx+q
constr_param.m(1)   =  0; % zero for standard case
constr_param.q(1)   = 15;

% Lower bound y<mx+q
constr_param.m(2)   =   0; % zero for standard case
constr_param.q(2)   =   0;



%% initial states
% Tractor
xt      =  0;                   % tractor inertial X position (m)
yt      =  0;                   % tractor inertial Y position (m)
psit    =    pi/2;              % tractor yaw angle (rad)
vt      =    4/3.6;             % tractor body x velocity (m/s) 

% Implement 
xi      =   0;                  % implement inertial X position (m)
yi      =   0;                  % implement inertial Y position (m)
psii    =   psit;               % implement yaw angle (rad)
vi      =   vt;                 % implement body x velocity (m/s)

z0      =   [xt;yt;psit;vt];

if strcmp(MODE,'01')
    % In case implement is attached, initial states are augmented and modified accordingly 
    
    z0=[xi+Li*cos(psii);yi+Li*sin(psii);psit;vt;xi;yi;psii;vi];
end



%% final state
% Tractor
xtf     =   xt + d;                                         % tractor inertial X position (m)
ytf     =   constr_param.m(2)*xtf + constr_param.q(2);      % tractor inertial Y position (m)
psitf   =   -pi/2;                                          % tractor yaw angle (rad)
vtf     =   4/3.6;                                          % tractor body x velocity (m/s) 

% Implement
xif     =    xi+d;                                          % implement inertial X position (m)
yif     =    constr_param.m(2)*xif + constr_param.q(2);     % implement inertial Y position (m)
psiif   =    psitf;                                         % implement yaw angle (rad)
vif     =    vt;                                            % implement body x velocity (m/s)

zf      =    [xtf;ytf;psitf;vtf];

if strcmp(MODE,'01')
     % In case implement is attached, final states are augmented and modified accordingly 

     zf=[xif+Li*cos(psiif); yif+Li*sin(psiif); psitf;vtf;xif;yif;psiif;vif];
end

constr_param.zf = zf;
%% Control problem parameters

Ns          =   75;                 % Simulation steps
Ts          =   0.25;               % initial guess for time step
Nu          =   2;                  % Down-sampling quantization for input variables

% Input saturation
vsat        =   15/3.6;             % Maximum tractor velocity (m/s)
asat        =   1;                  % Maximum tractor acceleration (m/s)
deltasat    =   30*pi/180;          % Maximum tractor steering angle (rad)                
delta_psi_sat = 75*pi/180;          % Maximum relative angle between tractor and implement (rad)

% Tolerances for the final state error
tol_f = ... 
    [0.05,0.05,5*pi/180,0.5/3.6,... % Tractor states tolerances
    0.05,0.05,5*pi/180,0.5/3.6]';   % Implement states tolerances

constr_param.vsat           =   vsat;
constr_param.delta_psi_sat  =   delta_psi_sat;
constr_param.tol_f          =   tol_f;

Np=ceil((Ns+1)/Nu);                 % Finite Horizon for the FHOCP
solution_flag=0;                    % Index to evaluate if the overall iteration output a valid solution or not

%% Linear Constraints
% Initialized vector accordingly to fmincon linear constraints notation 
% lb <= x <= ub 

lb       =       [-deltasat*ones(Np,1);
                 -asat*ones(Np,1)];

ub        =        [deltasat*ones(Np,1);
                   asat*ones(Np,1)];
 
%% Matlab fmincon options

options = optimoptions(@fmincon,...
    'Algorithm','interior-point',...
    'FiniteDifferenceType','central',...
    'MaxFunctionEvaluations',1e6, ...
    'MaxIterations',250,...
    'StepTolerance',1e-8,...
    'OptimalityTolerance', 10e-10,...
    'HessianApproximation', 'bfgs', ...
    'PlotFcn', {@plotfun_tractor_states},... 
    'Display','iter');

%% Run solver - Bulb trajectory 
tic ;

constr_param.c_vel = 0;                 % Force tractor to proceed forward, with positive velocities only          

U0     = [0.5*ones(8,1);                % Initial input sequence
           -0.5*ones(Np ...
           -8,1); 
           0.2*ones(ceil(Np/2),1);
           -0.2*ones(floor(Np/2),1);
           Ts;];                        % Include sampling time as optimization variable 

% Uncomment this line to evaluate initial condition trajectory
% [z01] = Tractor_traj(U0,z0,Nu,Ns,parameters,MODE);


[Ustar,fxstar,niter,exitflag,xsequence] = fmincon(@(U)cost_tractor_mincon(U,z0,Nu,Ns,parameters,constr_param,MODE)...
                                                    ,U0,[],[],[],[],lb,ub,...
                                                    @(U)constr_tractor_mincon(U,z0,Nu,Ns,parameters,constr_param,MODE),options);


%% Run solver - Fishtail trajectory 
% If no feasible solution has been found in the previous section, another
% optimization routine is initialized and computed, with different settings

if exitflag.constrviolation >options.ConstraintTolerance 

     constr_param.c_vel = 1;            % Include also negative velocities, to make the tractor exploit reverse gear          

    U0     = [-0.5*ones(12,1);          % Initial input sequence
                0.5*ones(12,1);
               -0.5*ones(Np-24,1); 
               0.2*ones(10,1);
               -0.7*ones(Np-27,1);
               0.5*ones(17,1);
               Ts;];                    % Include sampling time as optimization variable 

    
    % Uncomment this line to evaluate initial condition trajectory
    % [z02] = Tractor_traj(U0,z0,Nu,Ns,parameters,MODE);

    [Ustar,fxstar,niter,exitflag,xsequence] = fmincon(@(U)cost_tractor_mincon(U,z0,Nu,Ns,parameters,constr_param,MODE)...
                                                    ,U0,[],[],[],[],lb,ub,...
                                                    @(U)constr_tractor_mincon(U,z0,Nu,Ns,parameters,constr_param,MODE),options);

    if exitflag.constrviolation >options.ConstraintTolerance
        solution_flag=1;
        disp('No feasible solution has been found');
    end
end

opt_routine_time = toc;
results.opt_routine_time = opt_routine_time;

%% Compute and visualize final results
% If a valid solution it's found
if solution_flag==0
    % zstar is the evolution of the states applying the best input sequence U0 
    [zstar]         =   Tractor_traj(Ustar,z0,Nu,Ns,parameters,MODE);
     
    Ts_p            =   Ustar(end,1)*Nu;  
    Ts              =   Ustar(end,1); 
    
    results.Tend    =   Ts*Ns;
    results.Ts      =   Ts; 

    % Display Time variables
    fprintf('\n\n------------------------ Results -----------------------------\n\n');
    disp(['Final Time Tend: ', num2str(Ts*Ns), ' secondi']);
    disp(['Sampling Time Ts: ', num2str(Ts), ' secondi']);
    disp(['Optimization Routine Time: ', num2str(opt_routine_time), ' secondi']);
    
    % Export states and input for plotting
    plx             =   zstar(1,:)';
    ply             =   zstar(2,:)';

    ang             =   zstar(3,:)';
    vel             =   zstar(4,:)';
    
    delta           =   Ustar(1:Np,1);
    acc             =   Ustar(Np+1:end-1,1);
    
    asse            =   linspace(-5,10,2);

    close all 

    %% Save figures and plot results
    % Uncomment this to locally save images 
    % m_string    = replace(sprintf('m%.1f',constr_param.m(1)),'.','');
    % figname     = sprintf('___%s__q%.f',m_string,constr_param.q(1));
    % 
    % setup_name  = sprintf('\\%s_%s',MODE,figname);
    % mkdir([pwd '\Images' setup_name]);
    
    % X-Y trajectory plot
    fig1 = figure(1);
    fig1.OuterPosition(3:4) = [420 620]; 
    plot(plx,ply,'Color','b','DisplayName', 'Tractor', 'linewidth',1);hold on;
    plot(asse,constr_param.m(2)*asse + constr_param.q(2),"r",'DisplayName','Lower limit' ,'linewidth',1); hold on;
    plot(asse,constr_param.m(1)*asse + constr_param.q(1),"r",'DisplayName', 'Upper limit','linewidth',1); hold on;
    daspect([1 1 1]); axis([-5 10 -5 constr_param.q(1)+5]);
    xlabel('X [m]'); ylabel('Y [m]', 'Rotation', 90),grid on

    xLimits = get(gca, 'XLim');
    yLimits = get(gca, 'YLim');
    set(gca, 'XTick', xLimits(1):5:xLimits(2));
    set(gca, 'YTick', yLimits(1):5:yLimits(2))
    legend('show');

    % Input variables plot
    fig2 = figure(2);   
    fig2.OuterPosition(3:4) = [520 310]; 
    subplot(2,1,1);plot(0:Ts_p:(Np-1)*Ts_p,delta,'b','linewidth',1),xlabel('Time [s]'),ylabel('\delta_t  [rad]', 'Rotation', 90),grid on
    subplot(2,1,2);plot(0:Ts_p:(Np-1)*Ts_p,acc,'b','linewidth',1),xlabel('Time [s]'),ylabel('a_t [m/s^2]', 'Rotation', 90),grid on;
        
    if strcmp(MODE,'00')
        % Add final target point to the trajectory figure
        figure(fig1);
        plot(zf(1),zf(2),"xr",'MarkerSize', 10, 'LineWidth', 2,'DisplayName', 'Target point');
        legend('Location','northwest');
        
        % States plot 
        fig3 = figure(3);
        fig3.OuterPosition(3:4) = fig2.OuterPosition(3:4);
        subplot(2,1,1),plot(0:Ts:Ns*Ts,ang,'b','linewidth',1),xlabel('Time [s]'),ylabel('\psi_t [rad]', 'Rotation', 90),grid on
        
        subplot(2,1,2),plot(0:Ts:Ns*Ts,vel,'b','linewidth',1),xlabel('Time [s]'),ylabel('v_t [m/s]', 'Rotation', 90),grid on
    
    end
    
    if strcmp(MODE,'01')
        delta_psi=abs(zstar(3,:)-zstar(7,:));
        % Add final target point and implement trajectory to the figure
        figure(fig1);
        plot(zstar(5,:),zstar(6,:),'Color','g','DisplayName', 'Implement', 'linewidth',1); % hold on
        plot(zf(5),zf(6),"xr",'MarkerSize', 10, 'LineWidth', 2,'DisplayName', 'Target point');
        legend('Location','northwest');

        fig3 = figure(3);
        fig3.OuterPosition(3:4) = [720 310];
        subplot(2,2,1),plot(0:Ts:Ns*Ts,ang,'b','linewidth',1),xlabel('Time [s]'),ylabel('\psi_t [rad]', 'Rotation', 90),grid on
        
        subplot(2,2,3),plot(0:Ts:Ns*Ts,vel,'b','linewidth',1),xlabel('Time [s]'),ylabel('v_t [m/s]', 'Rotation', 90),grid on
         
        subplot(2,2,2),plot(0:Ts:Ns*Ts,zstar(7,:),'g','linewidth',1),xlabel('Time [s]'),ylabel('\psi_i [rad] ', 'Rotation', 90),grid on
        
        subplot(2,2,4),plot(0:Ts:Ns*Ts,zstar(8,:),'g','linewidth',1),xlabel('Time [s]'),ylabel('v_i [m/s]', 'Rotation', 90),grid on


    
        fig4 = figure(4);
       
        plot(0:Ts:(Ns)*Ts,delta_psi,'linewidth',1,'DisplayName', 'Relative angle'); hold on
        plot(0:Ts:(Ns)*Ts,delta_psi_sat+0*(0:Ts:(Ns)*Ts),'linewidth',1,'LineStyle','--','DisplayName', 'Maximum error')
        xlabel('Time [s]'); ylabel('Relative angle [rad]', 'Rotation', 90),grid on
        legend('show');

        % Uncomment this to locally save images 
        % saveas(fig4,[pwd '\Images' setup_name '/DeltaImplement.svg'])


    
    end
    
    % Uncomment this to locally save images 
    % saveas(fig1,[pwd '\Images' setup_name '\Trajectory.svg'])
    % saveas(fig2,[pwd '\Images' setup_name '\InputVariables.svg'])
    % saveas(fig3,[pwd '\Images' setup_name '\States.svg'])
    
    results.err_tractor_X   =   abs(plx(end,1)-zf(1));
    results.err_tractor_Y   =   abs(ply(end,1)-zf(2));
    results.err_tractor_psi =   abs(ang(end,1)-psitf);
    results.err_tractor_v   =   abs(vel(end,1)-vtf);
    if MODE=='01'
        results.err_implement_X =   abs(zstar(5,end)-zf(5));
        results.err_implement_Y =   abs(zstar(6,end)-zf(6));
        results.err_implement_psi = abs(zstar(7,end)-psiif);
        results.err_implement_v =   abs(zstar(8,end)-vif);
    end
    % Uncomment this to locally save results in a csv file 
    % T= struct2table(results);
    % writetable(T, [pwd '\Images' setup_name '\Results.csv']);
    
    
    fprintf('\n%-15s %-15s %-15s %-15s\n', 'Variable', 'Final Value', 'Target Value', 'Error');
    fprintf('-------------------------------------------------------------\n');
    fprintf('%-15s %-15.4f %-15.4f %-15.4f\n', 'xt', plx(end,1), zf(1), results.err_tractor_X);
    fprintf('%-15s %-15.4f %-15.4f %-15.4f\n', 'yt', ply(end,1), zf(2), results.err_tractor_Y);
    fprintf('%-15s %-15.4f %-15.4f %-15.4f\n', 'psit', ang(end,1), psitf, results.err_tractor_psi);
    fprintf('%-15s %-15.4f %-15.4f %-15.4f\n\n', 'vt', vel(end,1), vtf, results.err_tractor_v);
    
    if strcmp(MODE, '01')
        fprintf('%-15s %-15.4f %-15.4f %-15.4f\n', 'xi', zstar(5,end), zf(5), results.err_implement_X);
        fprintf('%-15s %-15.4f %-15.4f %-15.4f\n', 'yi', zstar(6,end), zf(6), results.err_implement_Y);
        fprintf('%-15s %-15.4f %-15.4f %-15.4f\n', 'psii', zstar(7,end), psiif, results.err_implement_psi);
        fprintf('%-15s %-15.4f %-15.4f %-15.4f\n', 'vi', zstar(8,end), vif, results.err_implement_v);
    end

else 
    clc
    close all
    disp('No feasible solution has been found')
end





