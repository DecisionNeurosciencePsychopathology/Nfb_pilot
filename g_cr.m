function [gx] = g_cr(Xt,Phi,u_t,inG)
% Identity observation mapping (partially observable)

beta = exp(Phi(1));
cs_next = u_t(3);
if cs_next || inG.noCS
[gx,dGdx,dGdP] = sigm(Xt(1),[],Phi);
else
[gx,dGdx,dGdP] = sigm(Xt(2),[],Phi);
end


% 
% if size(Phi,1) > 0
%     dG_dPhi = zeros(size(Phi,1),size(G,1));
% else
%     dG_dPhi = [];
% end
