function  [ fx] = f_pavlov( x_t,theta,u_t,in )
% function  [ fx,dfdx,dfdP ] = f_pavlov( x_t,P,u_t,in )
% IN:
% - x_t : Q-values (2*1)
% - P : learning rate (1*1)
% - u_t : previous action and feedback
% - in : []

alpha = 1./(1+exp(-theta(1)));          % learning rate
% epsilon =  1./(1+exp(-theta(2)));  % feedback decay
epsilon = 1;
gamma =  alpha;

cs = u_t(3);
us = u_t(4);
congruent = u_t(6);
feedback = u_t(4);

% feedback_salience = epsilon*(x_t(3));
if in.decay
feedback_salience = epsilon*x_t(3);
else
feedback_salience = 1;
end

if in.noCS
    fx(1) = x_t(1) + alpha*(feedback_salience*feedback-x_t(1));
else 
    if (cs && congruent) || (~cs && ~congruent) % actual infusion on previous trial
        fx(1) = x_t(1) + alpha*(feedback_salience*feedback-x_t(1));
        fx(2) = x_t(2)+ (alpha-gamma)*(feedback_salience*feedback-x_t(1));
    elseif (~cs && congruent) || (cs && ~congruent)  % no actual infusion on previous trial
        fx(2) = x_t(2) + alpha*(feedback_salience*feedback-x_t(2));
        fx(1) = x_t(1)+ (alpha-gamma)*(feedback_salience*feedback-x_t(2));
    end
end


if in.decay
fx(3) = (feedback_salience);  % this is the extent to which feedback has decayed 
end

%
% if cs
% %dfdx = [df1dx1 , df1dx2;
% %        df2dx1 , df2dx2]
%     dfdx = [1-alpha];
%     dfdP = [alpha*(1-alpha)*(us-x_t)];
%
%
% else
%     dfdx = 0;
%     dfdP = 0;
%
% end

% dfdx = dfdx';
% dfdP = dfdP';