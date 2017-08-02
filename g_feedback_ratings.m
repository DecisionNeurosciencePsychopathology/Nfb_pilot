function [gx] = g_feedback_ratings(Xt,Phi,u_t,inG)
% Identity observation mapping (partially observable)

beta = exp(Phi(1));  %% temperature
eta = 1./(1+exp(-Phi(2)));  %% expectancy sensitivity
epsilon = 1./(1+exp(-Phi(3)));  %% neurofeedback sensitivity
%lambda = Phi(4)./10;  %% bias
lambda = 0;
cs = u_t(3);
if cs || inG.noCS
[gx,~,~] = sigm((Xt(1)*eta + u_t(2)*epsilon + lambda),[],beta);
% exp_vigor = eta*Xt(1);
else
[gx,~,~] = sigm((Xt(2)*eta + u_t(2)*epsilon + lambda),[],beta);
% exp_vigor = eta*Xt(2);
end


% 
% if size(Phi,1) > 0
%     dG_dPhi = zeros(size(Phi,1),size(G,1));
% else
%     dG_dPhi = [];
% end
