function f = Tractor_cost_alternative(U,z0,zf,Np,Ns,parameters,Optimization_opt)

Q = Optimization_opt.Q;         
Qf = Optimization_opt.Qf;
Qdot = Optimization_opt.Qdot;
R = Optimization_opt.R;
vsat = Optimization_opt.vsat;
deltasat = Optimization_opt.deltasat;
asat = Optimization_opt.asat;
alpha = Optimization_opt.alpha;
beta = Optimization_opt.beta;
Ts_s=Optimization_opt.Ts_s;
Ts_p=Optimization_opt.Ts_p;
Tend=Optimization_opt.Tend;  

u_in        =   [U(1:Np,1)';
                U(Np+1:end,1)'];




%% Simulation 
zdot        =   zeros(4,1);
e           =   zeros(4,1);
ztemp       =   z0; 
f = 0;

     for ind=2:Ns+1
            u               =  u_in(:,ceil((ind-1)*Ts_s/Ts_p));
            zdot    =  tractor_model (ztemp,u,parameters);
            ztemp    =  ztemp+Ts_s*zdot;
            e        =  zf - ztemp;
    
    
    %% Computation of the cost function
    
    f   =   f + ...
            e'*Q*e + ...
            zdot'*Qdot*zdot + ...
            u'*R*u + ...
            alpha*exp(-beta*(0+ztemp(4,1))) + ...
            alpha*exp(-beta*(vsat-ztemp(4,1))) + ...
            alpha*exp(-beta*(deltasat+u(1))) + ...
            alpha*exp(-beta*(deltasat-u(1))) + ...
            alpha*exp(-beta*(asat+u(2))) +...
            alpha*exp(-beta*(asat-u(2)));
    
    
    end

f = f + e'*Qf*e;


end
