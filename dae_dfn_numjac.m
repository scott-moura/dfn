function ydot = dae_dfn_numjac(t,y,p)

x = y(1:size(p.f_x,1));
z = y(size(p.f_x,1)+1 : end);

[f,g] = dae_dfn(x,z,0,p);
ydot = [f; g];