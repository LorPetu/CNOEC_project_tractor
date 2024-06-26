function f = Tractor_cost_alternative_2(U,z0,zf,Ts,Np,parameters,Optimization_opt)

Q = Optimization_opt.Q;         
Qf = Optimization_opt.Qf;
Qdot = Optimization_opt.Qdot;
R = Optimization_opt.R;
vsat = Optimization_opt.vsat;
deltasat = Optimization_opt.deltasat;
asat = Optimization_opt.asat;
alpha = Optimization_opt.alpha;
beta = Optimization_opt.beta;

%% Build vector of inputs

t_in        =   [0:Ts:(Np-1)*Ts]';
u_in        =   [U(1:Np,1)';
                U(Np+1:end,1)'];

assignin('base','z0',z0);
assignin('base','t_in',t_in);
assignin('base','u_in',u_in);


%% Simulation 

Ts_under = Optimization_opt.Ts_under;

time_FFD    =   [0:Ts_under:(Np-1)*Ts];
Nblock      =   Ts/Ts_under;
Nsim_FFD    =   length(time_FFD);

z_sim       =   zeros(4,Nsim_FFD);
zdot        =   zeros(4,Nsim_FFD);
e           =   zeros(4,Nsim_FFD);
z_sim(1:4,1) =   z0;
f = 0;

    for ind=2:Nsim_FFD
    
        u               =  u_in(:,1+floor(time_FFD(ind)/Ts));
        zdot(:,ind)     =  Tractor_model (z_sim(:,ind-1),u,parameters);
        z_sim(:,ind)    =  z_sim(:,ind-1)+Ts/Nblock*zdot(:,ind);
        e(:,ind)        =  zf - z_sim(:,ind-1);
    
    
    %% Computation of the cost function
    
    f   =   f + ...
            e(:,1)'*Q*e(:,1) + ...
            zdot(:,1)'*Qdot*zdot(:,1) + ...
            u(:,1)'*R*u(:,1) + ...
            alpha*exp(-beta*(vsat+z_sim(4,ind))) + ...
            alpha*exp(-beta*(vsat-z_sim(4,ind))) + ...
            alpha*exp(-beta*(deltasat+u(1,1))) + ...
            alpha*exp(-beta*(deltasat-u(1,1))) + ...
            alpha*exp(-beta*(asat+u(2,1))) +...
            alpha*exp(-beta*(asat-u(2,1)));
    
    
    end

f = f + e(:,end)'*Qf*e(:,end);


end
