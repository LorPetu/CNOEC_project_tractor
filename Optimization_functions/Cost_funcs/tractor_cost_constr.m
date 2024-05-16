function v = tractor_cost_constr(U,z0,zf,parameters,Optimization_opt)
% Function that computes the trajectory of the tractor exiting from a row
% and the constraint

vsat = Optimization_opt.vsat;
deltasat = Optimization_opt.deltasat;
asat = Optimization_opt.asat;
Ts_s=Optimization_opt.Ts_s;
Ts_p=Optimization_opt.Ts_p;
Tend=Optimization_opt.Tend;  
Ns = Optimization_opt.Ns;
Np = Optimization_opt.Np;



u_in        =   [U(1:Np,1)';
                U(Np+1:end,1)'];

%% Run simulation with FFD
zdot        =   zeros(4,1);
e           =   zeros(4,1);
ztemp       =   z0; 
z_sim      =   zeros(4,Ns+1);
z_sim(:,1) =   z0;
f=0;
t=1;
for ind=2:Ns+1
    if abs(z_sim(1,ind-1) -zf(1))>0.025 || abs(z_sim(2,ind-1) -zf(2))>0.025 ||abs(z_sim(3,ind-1) -zf(3))>0.05 || abs(z_sim(4,ind-1) -zf(4))>0.05
        u               =  u_in(:,ceil((ind-1)*Ts_s/Ts_p));
        zdot               =   tractor_model(z_sim(:,ind-1),u,parameters);
        z_sim(:,ind)       =   z_sim(:,ind-1)+Ts_s*zdot;
        e = z_sim(1:2, ind)-z_sim(1:2, ind-1);
        
       f=f-1e-2*(z_sim(4,ind));
        
        t=t+1;
    else
        z_sim(:,ind)=zf;
        e = z_sim(1:2, ind)-z_sim(1:2, ind-1);
        
       
    end
     f= f+100*(e'*e);
end 

delta_delta=u_in(1,2:end)-u_in(1,1:end-1);
delta_acc=u_in(2,2:end)-u_in(2,1:end-1);
f=f+1e3*t;

%f=f+1e-6*(delta_acc*delta_acc');% 1e-5*(z_sim(4,:)*z_sim(4,:)');

   %così anche dopo che è arrivato tutto dovrebbe essere messo a zero           


%% Equality constraints g(x)
g = [z_sim(:,t)-zf];

%% Inequality constraints h(x)

h = [(z_sim(4,2:end)+vsat*ones(1,Ns))'; %not zero, otherwhise constraint would not be satisfied
    (-z_sim(4,2:end)+vsat*ones(1,Ns))'];
    %(-z_sim(2,2:end)+6.0*ones(1,Ns))'];
    %(z_sim(2,2:end)+0*ones(1,Ns))'];

%% Stack cost and constraints
v           =   [f;g;h];

end