function [z_sim] = stati_finali(U,z0,zf,Np,Ns,parameters,Optimization_opt)


vsat = Optimization_opt.vsat;
deltasat = Optimization_opt.deltasat;
asat = Optimization_opt.asat;
Ts_s=Optimization_opt.Ts_s;
Ts_p=Optimization_opt.Ts_p;
Tend=Optimization_opt.Tend; 
%% Build vector of inputs

u_in        =   [U(1:Np,1)';
                U(Np+1:end,1)'];

%% Simulation 
Ts_s=Optimization_opt.Ts_s;
Ts_p=Optimization_opt.Ts_p;
Tend=Optimization_opt.Tend; 

ztemp=z0;
z_sim      =   zeros(4,Ns);
z_sim(:,1) =   z0;


  for ind=2:Ns+1
    
        if abs(z_sim(1,ind-1) -zf(1))>0.2 || abs(z_sim(2,ind-1) -zf(2))>0.2 ||abs(z_sim(3,ind-1) -zf(3))>0.2 || abs(z_sim(4,ind-1) -zf(4))>0.2
        u               =  u_in(:,ceil((ind-1)*Ts_s/Ts_p));
        zdot               =   tractor_model(z_sim(:,ind-1),u,parameters);
        z_sim(:,ind)       =   z_sim(:,ind-1)+Ts_s*zdot;
        
        else
        z_sim(:,ind)=zf;
        end
  end
end