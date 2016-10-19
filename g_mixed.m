function [gx] = g_mixed(Xt,Phi,u_t,inG)
% Identity observation mapping (partially observable)

beta = exp(Phi(1));  %% temperature
eta = 1./(1+exp(-Phi(2)));  %% expectancy sensitivity
% epsilon = 1./(1+exp(-Phi(3)));  %% neurofeedback sensitivity
epsilon = 1 ; % fixed NF sensitivity
lambda = Phi(3)./10;  %% bias
% lambda = 0;
cs = u_t(3);
if cs || inG.noCS
[gx(1),~,~] = sigm((Xt(1)*eta + u_t(2)*epsilon + lambda),[],beta); %% prediction of feedback ratings as a function of 
...expectancy (Xt) weighted by expectancy sensitivity eta + actual neurofeedback presented (u_t(2)) + bias lambda
    
[gx(2),~,~] = sigm((Xt(1)),[],beta); %% prediction of infusion/no infusion expectancy ratings as a function of expectancy (Xt)

% exp_vigor = eta*Xt(1);
else
[gx(1),~,~] = sigm((Xt(2)*eta + u_t(2)*epsilon + lambda),[],beta);
[gx(2),~,~] = sigm((Xt(2)),[],beta);

% exp_vigor = eta*Xt(2);
end


% 
% if size(Phi,1) > 0
%     dG_dPhi = zeros(size(Phi,1),size(G,1));
% else
%     dG_dPhi = [];
% end
