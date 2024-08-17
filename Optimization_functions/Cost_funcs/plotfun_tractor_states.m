function stop = plotfun_tractor_states(U,optimValues,state)
persistent Ns Nu z0 parameters constr_param MODE% Retain these values throughout the optimization
global iterationCount;
stop = false;
switch state
    case "init"
        constr_param= evalin('base', 'constr_param');
        Nu = evalin('base', 'Nu');
        Ns = evalin('base', 'Ns');
        parameters = evalin('base', 'parameters');
        z0 = evalin('base', 'z0');
        MODE = evalin('base','MODE');

        zf=constr_param.zf;
        
        
      
    case "iter"
        %% Build vector of inputs
        zf=constr_param.zf;
        ztemp=z0;
        Np=ceil((Ns+1)/Nu);
        
        u_in        =   [U(1:Np,1)';
                        U(Np+1:2*Np,1)'];
        
        Ts=     U(end,1); 
        
        Tend=Ts*Ns;
        %% Simulate trajectory
        n_mode      = size(z0,1);
        
        zdot        =   zeros(n_mode,1);
        z_sim       =   zeros(n_mode,Ns+1);
        z_sim(:,1) =   z0;

        
        Tractor_model_used = str2func(['Tractor_',MODE, '_trail_model']);


        for ind=2:Ns+1
            u               =  u_in(:,ceil(ind/Nu));
            
            zdot               =   Tractor_model_used(z_sim(:,ind-1),u,parameters);
            z_sim(:,ind)       =   z_sim(:,ind-1)+Ts*zdot;
    
        end

        %% Plot
             
        time_s=linspace(0,Ts*Ns,Ns+1);
        time_p=linspace(0,Ts*Ns,Np);
        x_axis=linspace(-5,10,2);
       
        subplot(3,2,1:4)
        if strcmp(MODE,'00')
            plot(z_sim(1,:),z_sim(2,:),'b', ...
            x_axis,constr_param.m(2)*x_axis + constr_param.q(2),"r",...
            x_axis,constr_param.m(1)*x_axis + constr_param.q(1),"r",...
            zf(1),zf(2),"xr",'MarkerSize', 10),daspect([1,1,1]),grid on
            xlabel('x'), ylabel('y'),title(sprintf('Trajectory at k = %d\nTend= %f',optimValues.iteration, Tend))
        else
            plot(z_sim(1,:),z_sim(2,:),'b', ...
            z_sim(5,:),z_sim(6,:),'g',...
            x_axis,constr_param.m(2)*x_axis + constr_param.q(2),"r",...
            x_axis,constr_param.m(1)*x_axis + constr_param.q(1),"r",...
            zf(5),zf(6),"xr",'MarkerSize', 10),daspect([1,1,1]),grid on,%axis([-5 10 -2 12])

        end
        subplot(3,2,5),plot(0:Ts:(Np-1)*Ts,u_in(1,:),'b'),xlabel('Time (s)'),ylabel('delta'),grid on
        subplot(3,2,6),plot(0:Ts:(Np-1)*Ts,u_in(2,:),'b'),xlabel('Time (s)'),ylabel('acc'),grid on

    case "done"
      
end
end