function [RS]=simu1004(INS,opcode,tag);

% 随机数发生器初始化
rand('twister',sum(100*clock));

op=[];opname=[];
[op,opname]=getOP(opcode);

A=[];B=[];C=[];
if INS<11
    [A,B,C]=getInstance(INS,0,0);
else
    [A,B,C]=getInstance(INS,INS,10);
end

tt=0;gg=0;succ=0;
RS=[];
% 设定仿真次数
N=100;
for i=1:N
    [t,g,f,u,v,w]=permGA4(A,B,C,op);
    u
    v
    w
    %RS=[RS; t, g, f, 0, u, 0, v, 0, w];
    tt=tt+t;
    gg=gg+g;
    succ=succ+floor(f);
end
tt=tt/N;gg=gg/N;succ=succ/N;
disp('[avgTime  avgEvoGen  SuccRate]=');
disp([tt,gg,succ]);
%disp('Running records are:');
%disp(RS);
%save(tag,'RS','-ascii');
return





