function [l2_norm]=L2_norm_calculator(em,t)
dt=t(2)-t(1);
l2_norm=sqrt(sum(em.^2)*dt);
end