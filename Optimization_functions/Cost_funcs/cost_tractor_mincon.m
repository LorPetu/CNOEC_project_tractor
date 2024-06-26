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

lb_vel=constr_param.lb_vel;


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
e_acc=zeros(1,Np);
e_delta=zeros(1,Np);
for ind=2:Ns+1
    
    u               =  u_in(:,ceil(ind/Nu));
    zdot               =   Tractor_01_trail_model(z_sim(:,ind-1),u,parameters);
    z_sim(:,ind)       =   z_sim(:,ind-1)+Ts*zdot;

    %f=f+10*((z_sim(1:2, ind)-z_sim(1:2, ind-1))'*(z_sim(1:2, ind)-z_sim(1:2, ind-1)));
    
    e_acc(1,ind-1)=asat-abs(u(2));
    e_delta(1,ind-1)=deltasat-abs(u(1));
end 

delta_delta=u_in(1,2:end)-u_in(1,1:end-1);
delta_acc=u_in(2,2:end)-u_in(2,1:end-1);

Q=ones(1,8);
if lb_vel==0
    Q(3)=3;         %per velocità positiva q=3
    Q(7)=3;
elseif lb_vel==1
    Q=Q*1.3;
end
f = f +Q*s+80*Ts;   %  che va bene per v positva è =100
%+ 0.001*(delta_acc*delta_acc')+ 0.001*(delta_delta*delta_delta');

f=f*100;      %questo serve per scalare la funzione. Serve perchè fmincon non può settare i valori di linsearch e quindi con questo riusciamo a cambiarli (credo).
                %se è più alto la ricerca è più lenta ma più precisa es(50
                %o 100)sembrano funzionare bene
end