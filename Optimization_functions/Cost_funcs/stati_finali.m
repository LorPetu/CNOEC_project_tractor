function [z_sim] = stati_finali(U,z0,zf,Nu,Ns,parameters,Optimization_opt)


vsat = Optimization_opt.vsat;
deltasat = Optimization_opt.deltasat;
asat = Optimization_opt.asat;
%% Build vector of inputs

Np=ceil(Ns/Nu);

u_in        =   [U(1:Np,1)';
                U(Np+1:2*Np,1)'];

s =    U(2*Np+1:end-1,1);

Ts=     U(end,1);
%% Simulation 
ztemp=z0;
z_sim      =   zeros(4,Ns);
z_sim(:,1) =   z0;


 for ind=2:Ns+1
    if ceil(ind/Nu)<Np
        u               =  u_in(:,ceil(ind/Nu));
    end
    zdot               =   tractor_model(z_sim(:,ind-1),u,parameters);
    z_sim(:,ind)       =   z_sim(:,ind-1)+Ts*zdot;

  end
end