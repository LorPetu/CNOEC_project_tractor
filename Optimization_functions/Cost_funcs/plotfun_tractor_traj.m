function stop = plotfun_tractor_traj(U,optimValues,state)
persistent Ns Nu z0 parameters constr_param vsat asat deltasat % Retain these values throughout the optimization
stop = false;
switch state
    case "init"
        constr_param= evalin('base', 'constr_param');
        Nu = evalin('base', 'Nu');
        Ns = evalin('base', 'Ns');
        parameters = evalin('base', 'parameters');
        z0 = evalin('base', 'z0');

        zf=constr_param.zf;
        
      
    case "iter"
        %% Build vector of inputs
        
        zf=constr_param.zf;
        ztemp=z0;
        Np=ceil((Ns+1)/Nu);
        
        u_in        =   [U(1:Np,1)';
                        U(Np+1:2*Np,1)'];
        
        s =    U(2*Np+1:end-1,1);
        Ts=     U(end,1); 

        %% Simulate trajectory
        z_sim      =   zeros(4,Ns+1);
        z_sim(:,1) =   z0;

        for ind=2:Ns+1
            u               =  u_in(:,ceil(ind/Nu));
            zdot               =   tractor_model(z_sim(:,ind-1),u,parameters);
            z_sim(:,ind)       =   z_sim(:,ind-1)+Ts*zdot;
    
        end

        %% Plot
        % ang     = z_sim(3,:);
        % vel     = z_sim(4,:);
        % 
        % maxvsat = vsat*ones(Ns+1);
        % minvsat = -vsat*ones(Ns+1);
        % 
        % maxdeltasat = deltasat*ones(Np);
        % mindeltasat = -deltasat*ones(Np);
        % 
        % maxasat = asat*ones(Np);
        % minasat = -asat*ones(Np);
        % 
        % del     = u_in(1,:);
        % acc     = u_in(2,:);
        % 
        % 
        % plxf = zf(1,1);
        % plyf = zf(2,1);
        
        
        time_s=linspace(0,Ts*Ns,Ns+1);
        time_p=linspace(0,Ts*Ns,Np);
        x_axis=linspace(-5,10,2);
       
        plot(z_sim(1,:),z_sim(2,:),...
            x_axis,constr_param.m(2)*x_axis + constr_param.q(2),"green",...
            x_axis,constr_param.m(1)*x_axis + constr_param.q(1),"green",...
            zf(1),zf(2),"xr",'MarkerSize', 10),daspect([1,1,1]),grid on,axis([-5 10 -2 12])
        xlabel('x'), ylabel('y'),title(sprintf('Trajectory at k = %d',optimValues.iteration))
        



        % plot(time_s,ang),xlabel('Time (s)'),ylabel('psi'),
        % plot(time_s,vel*3.6,time_s,maxvsat*3.6,time_s,minvsat*3.6),xlabel('Time (s)'),ylabel('velocità [km/h]')
        % plot(time_p,del,time_p,maxdeltasat,time_p,mindeltasat),xlabel('Time (s)'),ylabel('delta ')
        % plot(time_p,acc,time_p,maxasat,time_p,minasat),xlabel('Time (s)'),ylabel('acc [m/s]');
        % 


    case "done"
        % Clean up plots
    % Some solvers also use case "interrupt"
end
end