function [tt,gg,succ]=simu1007(scale, max_gap)
% �������DDPʵ��

% �������������ʼ��
rand('twister',sum(100*clock));

op=[1];

tt=0;gg=0;succ=0;

% �趨�������
N=1000;

for i=1:N
    
    A=[];B=[];C=[];

    [A,B,C]=getInstance(100,scale,max_gap);
    
    [t,g,f,u,v,w]=permGA4(A,B,C,op);

    %disp([g,f]);
    %disp([g,f,0,u,0,v,0,w]);
    
    if floor(f)==1
        tt=tt+t;
        gg=gg+g;
        succ=succ+1;
    end

end

tt=tt/N;gg=gg/N;succ=succ/N;
disp('[avgTime  avgEvoGen  SuccRate]=');
disp([tt,gg,succ]);

return







