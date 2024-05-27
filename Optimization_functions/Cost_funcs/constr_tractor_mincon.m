function [h, g] = constr_tractor_mincon(U,z0,parameters,Optimization_opt, constr_param)
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

m_up= constr_param.m(1);
m_down= constr_param.m(2);
q_up = constr_param.q(1);
q_down = constr_param.q(2);

zf = constr_param.zf;




u_in        =   [U(1:Np,1)';
                U(Np+1:end,1)'];


%% Run simulation with FFD
zdot        =   zeros(4,1);
e           =   zeros(4,1);
ztemp       =   z0; 
z_sim      =   zeros(4,Ns+1);
z_sim(:,1) =   z0;

for ind=2:Ns+1
    u               =  u_in(:,ceil((ind-1)*Ts_s/Ts_p));
    zdot               =   tractor_model(z_sim(:,ind-1),u,parameters);
    z_sim(:,ind)       =   z_sim(:,ind-1)+Ts_s*zdot;
end 



%% Equality constraints g(x)
g = [z_sim(:,end)-zf];

%% Inequality constraints h(x)

%fmincon has the constraint structure c(x)<0, while myfmincon >0

h = [(-z_sim(4,2:end)+vsat*ones(1,Ns))'; %not zero, otherwhise constraint would not be satisfied
    (+z_sim(4,2:end)-vsat*ones(1,Ns))';
    % Parametrized constrained
    (z_sim(2,2:end)+m_up*z_sim(1,2:end)-q_up*ones(1,Ns))'; % -y + m*x + q < 0 y < m*x+q
    (-z_sim(2,2:end)+m_down*z_sim(1,2:end)+q_down*ones(1,Ns))']; % y - m*x - q > 0




end