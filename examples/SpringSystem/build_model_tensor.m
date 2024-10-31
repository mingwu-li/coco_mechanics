function [mass,damp,stiff,fnl,fext] = build_model_tensor()

n   = 2;
om1 = 2;
om2 = 4.5;
D1  = 0.01;
D2  = 0.2;
f1  = 0.02;
mass = eye(n,n);
damp = [2*D1*om1, 0;
     0, 2*D2*om2];
stiff = [om1^2 0;0 om2^2];
subs2 = [1 1 1
    1 2 2
    1 1 2
    2 2 2
    2 1 1
    2 1 2];
vals2 = [1.5*om1^2; 0.5*om1^2; om2^2; 1.5*om2^2; 0.5*om2^2; om1^2];
F2 = sptensor(subs2, vals2, [n,n,n]);
subs3 = [1 1 1 1
    1 1 2 2
    2 2 1 1
    2 2 2 2];
vals3 = 0.5*(om1^2+om2^2)*[1 1 1 1]';
F3  = sptensor(subs3, vals3, [n,n,n,n]);
fnl = {F2,F3};
fext = [f1;0];

end