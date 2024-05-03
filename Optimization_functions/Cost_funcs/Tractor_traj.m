function [z_sim] = Tractor_traj(U,z0,zf,Np,Ns,parameters,Optimization_opt)
%% Build vector of inputs

u_in        =   [   U(1:Np,1)';
                    U(Np+1:end,1)'];

%% Run simulation with FFD
Ts_s=Optimization_opt.Ts_s;
Ts_p=Optimization_opt.Ts_p;
Tend=Optimization_opt.Tend; 

ztemp=z0;
z_sim      =   zeros(4,Ns);
z_sim(:,1) =   z0;


  for ind=2:Ns
    
        u       =  u_in(:,ceil((ind-1)*Ts_s/Ts_p));
        zdot    =  tractor_model (ztemp,u,parameters);
        ztemp    =  ztemp+Ts_s*zdot;
        z_sim(:,ind)    =  ztemp;

  end

%% Post-processing the results

plx     = z_sim(1,:);
ply     = z_sim(2,:);
ang     = z_sim(3,:);
vel     = z_sim(4,:);

vsat = Optimization_opt.vsat;
maxvsat = vsat*ones(Ns);
minvsat = -vsat*ones(Ns);
deltasat = Optimization_opt.deltasat;
maxdeltasat = deltasat*ones(Np);
mindeltasat = -deltasat*ones(Np);
asat = Optimization_opt.asat;
maxasat = asat*ones(Np);
minasat = -asat*ones(Np);

del     = u_in(1,:);
acc     = u_in(2,:);


plxf = zf(1,1);
plyf = zf(2,1);

time_s=linspace(0,Ts_s,Ns);
time_p=linspace(0,Ts_p,Np);

figure(1);
subplot(5,2,1:6),plot(plx,ply,plxf,plyf,"xr"),daspect([1,1,1]),grid on,axis([-10 10 -10 10])
xlabel('x'), ylabel('y'),title('traiettoria')
subplot(5,2,7),plot(time_s,ang),xlabel('Time (s)'),ylabel('psi'),
subplot(5,2,8),plot(time_s,vel,time_s,maxvsat,time_s,minvsat),xlabel('Time (s)'),ylabel('velocità')
subplot(5,2,9);plot(time_p,del,time_p,maxdeltasat,time_p,mindeltasat),xlabel('Time (s)'),ylabel('delta')
subplot(5,2,10);plot(time_p,acc,time_p,maxasat,time_p,minasat),xlabel('Time (s)'),ylabel('acc');



end

