function v = tractor_cost_constr(U,z0,parameters,Optimization_opt, constr_param)
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
                U(Np+1:2*Np,1)'];

s =    U(2*Np+1:end,1);

%% Run simulation with FFD
zdot        =   zeros(4,1);
e_x          =   zeros(1,Ns);
e_y          =   zeros(1,Ns); 
ztemp       =   z0; 
z_sim      =   zeros(4,Ns+1);
z_sim(:,1) =   z0;
f=0;
peso=0;

for ind=2:Ns+1
    %peso=100*ind/(Ns+1);

    u               =  u_in(:,ceil((ind-1)*Ts_s/Ts_p));
    zdot               =   tractor_model(z_sim(:,ind-1),u,parameters);
    z_sim(:,ind)       =   z_sim(:,ind-1)+Ts_s*zdot;

    f=f+1e2*((z_sim(1:2, ind)-z_sim(1:2, ind-1))'*(z_sim(1:2, ind)-z_sim(1:2, ind-1)));
    
    e_x = (z_sim(1, ind)-zf(1));
    e_y = (z_sim(2, ind)-zf(2));

end 

delta_delta=u_in(1,2:end)-u_in(1,1:end-1);
delta_acc=u_in(2,2:end)-u_in(2,1:end-1);
delta_ex=e_x(2:end)-e_x(1:end-1);
delta_ey=e_y(2:end)-e_y(1:end-1);

f = f+ 1*(delta_acc*delta_acc')+ 1*(delta_delta*delta_delta')...
    +1e4*(s'*s);



    %+1e-1*(-(u_in(1,:)*u_in(1,:)'))...
    % +1e3*(e_x*e_x')...
    % +1e4*(e_y*e_y')...
    % 1e-4*(delta_acc*delta_acc')...
    %+1e-3*(delta_delta*delta_delta')...
    % +1*(delta_ex*delta_ex')...
    % +1*(delta_ey*delta_ey')...
    %+1e1*(-(u_in(1,:)*u_in(1,:)'))...
    %+1e8*(z_sim(:,end)-zf)'*(z_sim(:,end)-zf);

%% equality constraints 


g= [(abs(e_x(end))-s(1,1));                %uguaglianza con slask variables sulla fine
    (abs(e_y(end))-s(2,1));
    (abs(z_sim(3,end)-zf(3))-s(3,1));
    (abs(z_sim(4,end)-zf(4))-s(4,1))];
    
    
    % (z_sim(4,2:end)+ vsat*ones(1,Ns)+s(1:Ns,1)')';
    %  (-z_sim(4,2:end)+vsat*ones(1,Ns)+s(Ns+1:end,1)')'];


%% Inequality constraints h(x)

h = [(z_sim(4,2:end)+vsat*ones(1,Ns))'; %not zero, otherwhise constraint would not be satisfied
    (-z_sim(4,2:end)+vsat*ones(1,Ns))';
   
% Parametrized constrained
   (-z_sim(2,2:end)+m_up*z_sim(1,2:end)+(q_up+s(5,1))*ones(1,Ns))'; % -y + m*x + q > 0 y < m*x+q
   (z_sim(2,2:end)-m_down*z_sim(1,2:end)-q_down*ones(1,Ns))']; % y - m*x - q > 0


%% Stack cost and constraints

v           =   [f;g;h];

end