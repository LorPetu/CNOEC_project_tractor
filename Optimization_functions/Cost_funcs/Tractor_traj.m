function [z_sim] = Tractor_traj(U,z0,Nu,Ns,parameters,MODE)
    %% Build vector of inputs
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
        % FFD
        %zdot                =   Tractor_model_used(z_sim(:,ind-1),u,parameters);
        %z_sim(:,ind)        =   z_sim(:,ind-1)+Ts*zdot;
    
        % RK2
        z_prime_rk2         =   z_sim(:,ind-1)+Ts/2*Tractor_model_used(z_sim(:,ind-1),u,parameters);
        z_sim(:,ind)        =   z_sim(:,ind-1)+Ts*Tractor_model_used(z_prime_rk2,u,parameters);
    
    end


end

