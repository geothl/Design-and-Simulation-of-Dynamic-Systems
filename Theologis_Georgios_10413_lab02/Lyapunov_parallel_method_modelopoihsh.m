function dynamics = Lyapunov_parallel_method_modelopoihsh(dyn,u,a,b,Gamma,h)
dynamics(1)=-a*dyn(1)+b*u;
dynamics(2)=-dyn(3)*dyn(2)+dyn(4)*u;
dynamics(3)=-Gamma(1)*(dyn(1)+h-dyn(2))*dyn(2);
dynamics(4)=Gamma(2)*(dyn(1)+h-dyn(2))*u;
dynamics=transpose(dynamics);
end


% function dynamics = Lyapunov_parallel_method_modelopoihsh(dyn,u,a,b,Gamma,h)
% dynamics(1)=-a*dyn(1)+b*u;
% dynamics(2)=-dyn(3)*dyn(2)+dyn(4)*u;
% dynamics(3)=-Gamma(1)*(dyn(1)-dyn(2))*dyn(2);
% dynamics(4)=Gamma(2)*(dyn(1)-dyn(2))*u;
% dynamics=transpose(dynamics);
% end