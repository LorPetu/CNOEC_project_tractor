function f = cost_tractor_mincon(U,z0,parameters,Optimization_opt, constr_param)
% Function that computes the trajectory of the tractor exiting from a row
% and the constraint

vsat = Optimization_opt.vsat;
deltasat = Optimization_opt.deltasat;
asat = Optimization_opt.asat;

Ns = Optimization_opt.Ns;
Nu=Optimization_opt.Nu;

m_up= constr_param.m(1);
m_down= constr_param.m(2);
q_up = constr_param.q(1);
q_down = constr_param.q(2);

zf = constr_param.zf;




Np=ceil((Ns+1)/Nu);

u_in        =   [U(1:Np,1)';
                U(Np+1:2*Np,1)'];

s =    U(2*Np+1:end-1,1);

Ts=     U(end,1);

%% Run simulation with FFD

zdot        =   zeros(8,1);
z_sim      =   zeros(8,Ns+1);
z_sim(:,1) =   z0;
f=0;

for ind=2:Ns+1
    
    u               =  u_in(:,ceil(ind/Nu));
    zdot               =   Tractor_01_trail_model(z_sim(:,ind-1),u,parameters);
    z_sim(:,ind)       =   z_sim(:,ind-1)+Ts*zdot;

    %f=f+10*((z_sim(1:2, ind)-z_sim(1:2, ind-1))'*(z_sim(1:2, ind)-z_sim(1:2, ind-1)));

end 

delta_delta=u_in(1,2:end)-u_in(1,1:end-1);
delta_acc=u_in(2,2:end)-u_in(2,1:end-1);


f = f + 1e3*(s'*s)+1e4*Ts;
    %+ 1*(delta_acc*delta_acc')+ 1*(delta_delta*delta_delta');

end