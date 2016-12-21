function [gx] = g_mixed(Xt,Phi,u_t,inG)
% Identity observation mapping (partially observable)

beta1 = exp(Phi(1));  %% inverse temperature
beta2 = exp(Phi(2));  %% inverse temperature

% epsilon = 1./(1+exp(-Phi(3)));  %% neurofeedback sensitivity
% epsilon = 1 ; % fixed NF sensitivity
if inG.biases       %% allow for subject-level biases in feedback and expectancy ratings
    lambda1 = Phi(3)./10;  %% feedback bias
    lambda2 = Phi(4)./10;  %% expectancy bias
else
    lambda1 = 0;  %% feedback bias
    lambda2 = 0;  %% expectancy bias
end

% add trial-invariable infusion expectancy parameter
if inG.infusion_expectancy
    epsilon = Phi(3)./10;
else
    epsilon = 0;
end


% lambda = 0;

cs_current = u_t(1); % infusion cue on current trial
cs_previous = u_t(3); % infusion cue on previous trial
congruent_current = u_t(5); % current trial congruent
congruent_previous = u_t(6); % previous trial congruent

if inG.decay
if (cs_current && congruent_current) || (~cs_current && ~congruent_current) % trials where they received infusions
[gx(1),~,~] = sigm((Xt(1) + u_t(2)*(Xt(3)) + lambda1),[],beta1); %% prediction of feedback ratings as a function of 
... infusion expectancy (Xt) weighted by expectancy sensitivity eta + actual neurofeedback presented (u_t(2)) + bias lambda
    

% exp_vigor = eta*Xt(1);
elseif (~cs_current && congruent_current) || (cs_current && ~congruent_current) % trials with no infusions
[gx(1),~,~] = sigm((Xt(2) + u_t(2)*(Xt(3)) + lambda1),[],beta1);

% exp_vigor = eta*Xt(2);
end
else
if (cs_current && congruent_current) || (~cs_current && ~congruent_current) % trials where they received infusions
[gx(1),~,~] = sigm((Xt(1) + u_t(2) + lambda1 + epsilon),[],beta1); %% prediction of feedback ratings as a function of 
... infusion expectancy (Xt) weighted by expectancy sensitivity eta + actual neurofeedback presented (u_t(2)) + bias lambda
    

% exp_vigor = eta*Xt(1);
elseif (~cs_current && congruent_current) || (cs_current && ~congruent_current) % trials with no infusions
[gx(1),~,~] = sigm((Xt(2) + u_t(2) + lambda1),[],beta1);

% exp_vigor = eta*Xt(2);
end
end
    
    
if cs_current
[gx(2),~,~] = sigm((Xt(1) + lambda2 + epsilon),[],beta2); %% prediction of infusion/no infusion expectancy ratings as a function of expectancy (Xt)
else
[gx(2),~,~] = sigm((Xt(2) + lambda2),[],beta2);
end


% 
% if size(Phi,1) > 0
%     dG_dPhi = zeros(size(Phi,1),size(G,1));
% else
%     dG_dPhi = [];
% end
