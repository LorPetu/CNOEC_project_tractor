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
    
        u       =  u_in(:,ceil((ind-1)*Ts_s/Ts_p));
        zdot    =  tractor_model (ztemp,u,parameters);
        ztemp    =  ztemp+Ts_s*zdot;
        z_sim(:,ind)    =  ztemp;

  end
end