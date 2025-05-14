function dynamics = Lyapunov_seriesparallel_method_modelopoihsh(dyn,u,a,b,Gamma,h,thm)
dynamics(1)=-a*dyn(1)+b*u;
dynamics(2)=-dyn(3)*(dyn(1)+h)+dyn(4)*u+thm*(dyn(1)+h-dyn(2));
dynamics(3)=-Gamma(1)*(dyn(1)+h-dyn(2))*(dyn(1)+h);
dynamics(4)=Gamma(2)*(dyn(1)+h-dyn(2))*u;
dynamics=transpose(dynamics);
end