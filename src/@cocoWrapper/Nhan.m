function y = Nhan(obj,u,v)
% NHAN This function returns internal force N(u,v) in equation of motion as
% follows M\ddot{u}+N(u,\dot{u})=F(t,p). This function will be used in
% forward toolbox for finding periodic orbits
%
% Y = NHAN(OBJ,U,V)
%
% obj: coco object with dynamical system included
% u:   displacement
% v:   velocity
% 
% See also: DNDU, DNDV

y = obj.system.C*v+obj.system.K*u+obj.system.compute_fnl(u,v);

% y = obj.system.C*v+obj.system.K*u;
% another version for compute_dfnldx
% fnl = obj.multiFnl; % fnl.coeffs and fnl.ind
% numNonlinearTerms = size(fnl.coeffs,2);
% for i=1:numNonlinearTerms
%     coeff = fnl.coeffs(:,i);
%     ind = fnl.ind(i,:);
%     % find nonzero exponents
%     expind = find(ind);
%     s = 1;
%     for j=1:numel(expind)
%         s = s.*u(expind(j),:).^ind(expind(j));
%     end
%     y = y+coeff*s;
% end

end