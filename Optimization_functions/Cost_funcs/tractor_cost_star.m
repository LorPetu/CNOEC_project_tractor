function [f] = tractor_cost_star(U,z0,zf,Ts,Np,parameters,Optimization_opt)
[F] = tractor_cost_GN_grad(U,z0,zf,Ts,Np,parameters,Optimization_opt);


f=F'*F;
