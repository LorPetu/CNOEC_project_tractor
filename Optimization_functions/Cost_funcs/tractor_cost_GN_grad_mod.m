function F = tractor_cost_GN_grad_mod(U,z0,zf,Np,Ns,parameters,Optimization_opt)

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

%inizializzazione funzione di costo gauss- newton

F           =   [zeros(4*Ns,1);
                zeros(4*Ns,1);
                zeros(2*Ns,1);
                zeros(Ns,1);        %alpha*exp(-beta*(0+ztemp(4,1)));
                zeros(Ns,1);        %alpha*exp(-beta*(vsat+ztemp(4,1)));
                zeros(Ns,1);        %alpha*exp(-beta*(deltasat+U(1,1)));
                zeros(Ns,1);         %alpha*exp(-beta*(deltasat-U(1,1)));
                zeros(Ns,1);         %alpha*exp(-beta*(asat+U(2,1)));
                zeros(Ns,1);          %alpha*exp(-beta*(asat+U(2,1)));
                zeros(Ns,1);
                zeros(Ns,1);
                zeros(4,1);
                0];

t=0;

    for ind=2:Ns+1

        if abs(ztemp(1)-zf(1))>0.1 || abs( ztemp(2)-zf(2))>0.1
            u               =  u_in(:,ceil((ind-1)*Ts_s/Ts_p));
            zdot    =  tractor_model (ztemp,u,parameters);
            ztemp    =  ztemp+Ts_s*zdot;
            e        =  zf - ztemp;
       % Update the cost function 
        F((ind-2)*4+1:(ind-1)*4,1)          =   Q*[e(:,1)];
        F((ind-2)*4+4*Ns+1:(ind-1)*4+4*Ns,1)          =   Qdot*zdot;
        F((ind-2)*2+8*Ns+1:(ind-1)*2+8*Ns,1)              =   R*u;
        F(8*Ns+2*Ns+(ind-1),1)                    =   alpha*exp(-beta*(0+ztemp(4,1)));
        F(8*Ns+2*Ns+Ns+(ind-1),1)                    =   alpha*exp(-beta*(vsat-ztemp(4,1)));
        F(8*Ns+2*Ns+2*Ns+(ind-1),1)                    =   alpha*exp(-beta*(deltasat+u(1)));
        F(8*Ns+2*Ns+3*Ns+(ind-1),1)                    =   alpha*exp(-beta*(deltasat-u(1)));
        F(8*Ns+2*Ns+4*Ns+(ind-1),1)                    =   alpha*exp(-beta*(asat+u(2)));
        F(8*Ns+2*Ns+5*Ns+(ind-1),1)                    =   alpha*exp(-beta*(asat-u(2)));
        F(8*Ns+2*Ns+6*Ns+(ind-1),1)                    =   0;%alpha*exp(-beta*(0+ztemp(1,1)));
        F(8*Ns+2*Ns+7*Ns+(ind-1),1)                    =   0;%alpha*exp(-beta*(6-ztemp(2,1)));
        t=t+1;
        end
    end 
% Update the terminal cost and its Jacobian
F(8*Ns+2*Ns+8*Ns+1:8*Ns+2*Ns+8*Ns+4,1)                =  Qf*e;

if t<200
F(8*Ns+2*Ns+8*Ns+1:8*Ns+2*Ns+8*Ns+4+1,1)=1e-2*t;
end
disp(ztemp)

disp( t)
end


