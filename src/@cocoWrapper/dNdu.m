function y = dNdu(obj,u,v)
% DNDU This function returns the derivative of internal force N(u,v) with
% respect to its first argument. EOM is given by M\ddot{u}+N(u,\dot{u})=F(t,p). 
% This function will be used in forward toolbox for finding periodic orbits
%
% Y = DNDU(OBJ,U,V)
%
% obj: coco object with dynamical system included
% u:   displacement
% v:   velocity
% 
% See also: NHAN, DNDV
% tic
y = obj.system.K+obj.system.compute_dfnldx(u,v);
% toc
% tic
% another version for compute_dfnldx - incorrect in the case of zero s
% y = obj.system.K;
% fnl = obj.multiFnl;
% numNonlinearTerms = size(fnl.coeffs,2);
% for i=1:numNonlinearTerms
%     if i==184
%        disp('fsf'); 
%     end
%     coeff = fnl.coeffs(:,i);
%     ind = fnl.ind(i,:);
%     % find nonzero exponents
%     expind = find(ind);
%     numind = numel(expind);
%     s = 1;
%     for j=1:numind
%         s = s*u(expind(j))^ind(expind(j));
%     end
%     for j=1:numind
%         ij = ind(expind(j));
%         uj = u(expind(j));
%         if uj==0
%             ijoveruj = 0^(ij-1);
%         else
%             ijoveruj = ij/uj;
%         end
%         jac = ijoveruj*s;
%         y(:,expind(j)) = y(:,expind(j))+coeff*jac;
%     end   
% end
% toc

% Jx = coco_ezDFDX('f(x,p)', @obj.Nhan, u, v);

end