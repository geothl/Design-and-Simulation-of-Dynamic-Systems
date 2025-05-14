function dynamics = Lyapunov_seriesparallel_method_modelopoihsh_nonlinear(dyn,u,th1,th2,Gamma,thm,f)
dynamics(1)=-th1*f(dyn(1))+th2*u;
dynamics(2)=dyn(3)*f(dyn(1))+dyn(4)*u+thm*(dyn(1)-dyn(2));
dynamics(3)=-Gamma(1)*f(dyn(1))*(dyn(1)-dyn(2));
dynamics(4)=Gamma(2)*u*(dyn(1)-dyn(2));
dynamics=transpose(dynamics);

end