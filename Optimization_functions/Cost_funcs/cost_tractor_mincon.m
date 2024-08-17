function f = cost_tractor_mincon(U,z0,parameters,Optimization_opt, constr_param, MODE)
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

Ts=     U(end,1);

%% Run simulation with FFD

n_mode      = size(z0,1);

zdot        =   zeros(n_mode,1);
z_sim       =   zeros(n_mode,Ns+1);
z_sim(:,1)  =   z0;
f           =   0;
f1          =   0;
f2          =   0;
p           =   ones(size(zf));

p(3)        = 5;
%p(7)        = p(3);
Q           =   zeros(Ns,1);
Q(end,1)    =   1;
Tractor_model_used = str2func(['Tractor_',MODE, '_trail_model']);

for ind=2:Ns+1
    
    u               =  u_in(:,ceil(ind/Nu));
    zdot               =   Tractor_model_used(z_sim(:,ind-1),u,parameters);
    z_sim(:,ind)       =   z_sim(:,ind-1)+Ts*zdot;

    f1=f1+Q(ind-1,1)*p'*abs(z_sim(:,ind)-zf);
end 

f2 = Ns*Ts;

gamma =1;
f = gamma*f1 + (1-gamma)*f2; 
%disp(["f1 = ", num2str(f1),"f2= ",num2str(f2)]);

f=f*50;      %questo serve per scalare la funzione. Serve perchè fmincon non può settare i valori di linsearch e quindi con questo riusciamo a cambiarli (credo).
                %se è più alto la ricerca è più lenta ma più precisa es(50
                %o 100)sembrano funzionare bene
end