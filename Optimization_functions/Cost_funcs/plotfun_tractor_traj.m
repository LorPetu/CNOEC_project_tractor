function stop = plotfun_tractor_traj(U,optimValues,state)
persistent Ns Nu z0 parameters constr_param % Retain these values throughout the optimization
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
        
        Tend=Ts*Ns;
        %% Simulate trajectory
        zdot=zeros(8,1);
        z_sim      =   zeros(8,Ns+1);
        z_sim(:,1) =   z0;

        for ind=2:Ns+1
            u               =  u_in(:,ceil(ind/Nu));
            zdot               =   Tractor_01_trail_model(z_sim(:,ind-1),u,parameters);
            z_sim(:,ind)       =   z_sim(:,ind-1)+Ts*zdot;
    
        end

        %% Plot
              
        time_s=linspace(0,Ts*Ns,Ns+1);
        time_p=linspace(0,Ts*Ns,Np);
        x_axis=linspace(-5,10,2);
       
        plot(z_sim(1,:),z_sim(2,:),'blue', ...
            z_sim(5,:),z_sim(6,:),'black',...
            x_axis,constr_param.m(2)*x_axis + constr_param.q(2),"green",...
            x_axis,constr_param.m(1)*x_axis + constr_param.q(1),"green",...
            zf(5),zf(6),"xr",'MarkerSize', 10),daspect([1,1,1]),grid on,%axis([-5 10 -2 12])
        xlabel('x'), ylabel('y'),title(sprintf('Trajectory at k = %d\nTend= %f',optimValues.iteration, Tend))
        
    case "done"
      
end
end