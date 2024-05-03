function F = tractor_cost_GN_grad(U,z0,zf,Ts,Np,parameters,Optimization_opt)

Q = Optimization_opt.Q;         
Qf = Optimization_opt.Qf;
Qdot = Optimization_opt.Qdot;
R = Optimization_opt.R;
vsat = Optimization_opt.vsat;
deltasat = Optimization_opt.deltasat;
asat = Optimization_opt.asat;
alpha = Optimization_opt.alpha;
beta = Optimization_opt.beta;
%Ts_s=Optimization_opt.Ts_s;
%Ts_p=Optimization_opt.Ts_p;
%Tend=Optimization_opt.Tend;  

t_in        =   [0:Ts:(Np-1)*Ts]';
u_in        =   [U(1:Np,1)';
                U(Np+1:end,1)'];
%% Simulation 

Ts_under = Optimization_opt.Ts_under;

time_FFD    =   [0:Ts_under:(Np-1)*Ts];
Nblock      =   Ts/Ts_under;
Nsim_FFD    =   length(time_FFD);

z_sim       =   zeros(4,Nsim_FFD+1);
z_sim(:,1)  =   z0;
zdot        =   zeros(4,1);
e           =   zeros(4,1);
ztemp       =   z0; 
%inizializzazione funzione di costo gauss- newton


F           =   [zeros(4*Nsim_FFD,1);
                zeros(4*Nsim_FFD,1);
                zeros(2*Nsim_FFD,1);
                zeros(Nsim_FFD,1);        %alpha*exp(-beta*(0+ztemp(4,1)));
                zeros(Nsim_FFD,1);        %alpha*exp(-beta*(vsat+ztemp(4,1)));
                zeros(Nsim_FFD,1);        %alpha*exp(-beta*(deltasat+U(1,1)));
                zeros(Nsim_FFD,1);         %alpha*exp(-beta*(deltasat-U(1,1)));
                zeros(Nsim_FFD,1);         %alpha*exp(-beta*(asat+U(2,1)));
                zeros(Nsim_FFD,1);          %alpha*exp(-beta*(asat+U(2,1)));
                zeros(4,1);];



    for ind=2:Nsim_FFD
            u               =  u_in(:,1+floor(time_FFD(ind)/Ts));
            zdot    =  tractor_model (ztemp,u,parameters);
            ztemp    =  ztemp+Ts/Nblock*zdot;
            e        =  zf - ztemp;
            zsim(:,ind) = ztemp;
        
       % Update the cost function 
        F((ind-2)*4+1:(ind-1)*4,1)          =   Q*[e(:,1)];
        F((ind-2)*4+4*Nsim_FFD+1:(ind-1)*4+4*Nsim_FFD,1)          =   Qdot*zdot;
        F((ind-2)*2+8*Nsim_FFD+1:(ind-1)*2+8*Nsim_FFD,1)              =   R*u;
        F(8*Nsim_FFD+2*Nsim_FFD+(ind-1),1)                    =   alpha*exp(-beta*(0+ztemp(4,1)));
        F(8*Nsim_FFD+2*Nsim_FFD+Nsim_FFD+(ind-1),1)                    =   alpha*exp(-beta*(vsat-ztemp(4,1)));
        F(8*Nsim_FFD+2*Nsim_FFD+2*Nsim_FFD+(ind-1),1)                    =   alpha*exp(-beta*(deltasat+u(1)));
        F(8*Nsim_FFD+2*Nsim_FFD+3*Nsim_FFD+(ind-1),1)                    =   alpha*exp(-beta*(deltasat-u(1)));
        F(8*Nsim_FFD+2*Nsim_FFD+4*Nsim_FFD+(ind-1),1)                    =   alpha*exp(-beta*(asat+u(2)));
        F(8*Nsim_FFD+2*Nsim_FFD+5*Nsim_FFD+(ind-1),1)                    =   alpha*exp(-beta*(asat-u(2)));
   
    end 

% Update the terminal cost and its Jacobian
F(8*Nsim_FFD+2*Nsim_FFD+6*Nsim_FFD+1:8*Nsim_FFD+2*Nsim_FFD+6*Nsim_FFD+4,1)                =   Qf*e;

end


