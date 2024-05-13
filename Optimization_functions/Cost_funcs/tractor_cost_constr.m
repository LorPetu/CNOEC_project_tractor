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
for ind=2:Ns+1
    u               =  u_in(:,ceil((ind-1)*Ts_s/Ts_p));
    zdot               =   tractor_model(z_sim(:,ind-1),u,parameters);
    z_sim(:,ind)       =   z_sim(:,ind-1)+Ts_s*zdot;
    e = z_sim(1:2, ind)-z_sim(1:2, ind-1);

    f= f +1e2*(e'*e); 
end 

% delta_delta=u_in(1,2:end)-u_in(1,1:end-1);
delta_acc=u_in(2,2:end)-u_in(2,1:end-1);
% 
f=f+1*(delta_acc*delta_acc');%(delta_delta*delta_delta')+1e3*(delta_acc*delta_acc');


%% Equality constraints g(x)
g = [z_sim(:,end)-zf];

%% Inequality constraints h(x)

h = [(z_sim(4,2:end)+vsat*ones(1,Ns))'; %not zero, otherwhise constraint would not be satisfied
    (-z_sim(4,2:end)+vsat*ones(1,Ns))';
    (-z_sim(2,2:end)+6.0*ones(1,Ns))';
    (z_sim(2,2:end)+0*ones(1,Ns))'];

%% Stack cost and constraints
v           =   [f;g;h];

end