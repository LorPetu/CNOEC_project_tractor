function [f,zsim] = modified_tractor_cost_function(u, theta, z0, zf, Ts, Q, Qf, Qdot, R, alpha_barrier, beta_barrier)


[F,zsim] = modified_tractor_cost_function_F(u, theta, z0, zf, Ts, Q, Qf, Qdot, R, alpha_barrier, beta_barrier);

f=F'*F;





