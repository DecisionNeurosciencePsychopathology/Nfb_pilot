function  [ fx] = f_pavlov( x_t,theta,u_t,in )
% function  [ fx,dfdx,dfdP ] = f_pavlov( x_t,P,u_t,in )
% IN:
% - x_t : Q-values (2*1)
% - P : learning rate (1*1)
% - u_t : previous action and feedback
% - in : []

alpha = 1./(1+exp(-theta(1)));
% gamma =  1./(1+exp(-theta(2)));
gamma =  alpha;

cs = u_t(3);
us = u_t(4);
if in.noCS
    fx(1) = x_t(1) + alpha*(us-x_t(1));
else
    if cs
        fx(1) = x_t(1) + alpha*(us-x_t(1));
        fx(2) = x_t(2)+ (alpha-gamma)*(us-x_t(1));
    else
        fx(2) = x_t(2) + alpha*(us-x_t(2));
        fx(1) = x_t(1)+ (alpha-gamma)*(us-x_t(2));
    end
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