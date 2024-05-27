function stop = plotfun_tractor_traj(U,optimValues,state)
persistent Ts_s Ts_p Tend Ns Np z0 parameters vsat asat deltasat % Retain these values throughout the optimization
stop = false;
switch state
    case "init"
        %% Run simulation with FFD
        Ts_p=0.3; %Optimization_opt.Ts_p;
        Ts_s=0.2; %Optimization_opt.Ts_s;

        Tend=14; %Optimization_opt.Tend; 

        Ns          =   Tend/Ts_s;                    % Simulation steps
        Np          =   ceil(Tend/Ts_p);                  % Prediction steps
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
        
        %%
        vsat        =   20/3.6;                     % Input saturation
        asat        =   1;                      % Cart position limits
        deltasat    =   30*pi/180;
       
        
    case "iter"
        %% Build vector of inputs
        ztemp=z0;
        u_in        =   [   U(1:Np,1)';
                            U(Np+1:end,1)'];

        %% Simulate trajectory
        z_sim      =   zeros(4,Ns);
        z_sim(:,1) =   z0;

        for ind=2:Ns+1
        
            u       =  u_in(:,ceil((ind-1)*Ts_s/Ts_p));
            zdot    =  tractor_model (ztemp,u,parameters);
            ztemp    =  ztemp+Ts_s*zdot;
            z_sim(:,ind)    =  ztemp;
    
        end

        %% Plot
        ang     = z_sim(3,:);
        vel     = z_sim(4,:);

        maxvsat = vsat*ones(Ns+1);
        minvsat = -vsat*ones(Ns+1);
     
        maxdeltasat = deltasat*ones(Np);
        mindeltasat = -deltasat*ones(Np);
       
        maxasat = asat*ones(Np);
        minasat = -asat*ones(Np);

        del     = u_in(1,:);
        acc     = u_in(2,:);
        
        % 
        % plxf = zf(1,1);
        % plyf = zf(2,1);
        
        time_s=linspace(0,Tend,Ns+1);
        time_p=linspace(0,Tend,Np);
        plot(z_sim(1,:),z_sim(2,:)),daspect([1,1,1]),grid on,axis([-5 10 -2 12])
        xlabel('x'), ylabel('y'),title(sprintf('Trajectory at k = %d',optimValues.iteration))
        
        % plot(time_s,ang),xlabel('Time (s)'),ylabel('psi'),
        % plot(time_s,vel*3.6,time_s,maxvsat*3.6,time_s,minvsat*3.6),xlabel('Time (s)'),ylabel('velocità [km/h]')
        % plot(time_p,del,time_p,maxdeltasat,time_p,mindeltasat),xlabel('Time (s)'),ylabel('delta ')
        % plot(time_p,acc,time_p,maxasat,time_p,minasat),xlabel('Time (s)'),ylabel('acc [m/s]');



    case "done"
        % Clean up plots
    % Some solvers also use case "interrupt"
end
end