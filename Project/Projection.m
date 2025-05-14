function [para_ek_dot]= Projection(a_ek,b_ek,para_ek_dot,Gamma)

grad_g=0;


% if((a_ek>0)&&(b_ek>0.5)&&(b_ek<1.5))

% else

if((a_ek<=0)&&(b_ek>0.5)&&(b_ek<1.5))
 grad_g=[-1;0];   
elseif((a_ek>0)&&(b_ek<=0.5))
 grad_g=[0;-1];   
elseif((a_ek>0)&&(b_ek>=1.5))
 grad_g=[0;1];   
elseif((a_ek<=0)&&(b_ek>=1.5))
 grad_g=[-1/sqrt(2);1/sqrt(2)];   
elseif((a_ek<=0)&&(b_ek>=0.5))
 grad_g=[-1/sqrt(2);-1/sqrt(2)];   
 % else
 % para_ek_dot= Gamma*para_ek_dot;
end

% if(transpose(para_ek_dot)*grad_g>0)
% 
%  para_ek_dot= Gamma*para_ek_dot-Gamma*grad_g*transpose(grad_g)/(transpose(grad_g)*Gamma*grad_g)*Gamma*para_ek_dot;
% end


if((a_ek>0)&&(b_ek>0.5)&&(b_ek<1.5))
para_ek_dot= Gamma*para_ek_dot;
%includes area near border as border because of the nature of ode
elseif(((((a_ek<=0)&&(b_ek>0.5)&&(b_ek<1.5))||((a_ek>0)&&(b_ek<=0.5))||((a_ek>0)&&(b_ek>=1.5))))&&(transpose(para_ek_dot)*grad_g<=0))
para_ek_dot= Gamma*para_ek_dot;

else
para_ek_dot= Gamma*para_ek_dot-Gamma*grad_g*transpose(grad_g)/(transpose(grad_g)*Gamma*grad_g)*Gamma*para_ek_dot;
end

end