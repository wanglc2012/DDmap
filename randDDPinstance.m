function [A,B,C]=randDDPinstance(scale,max_gap)
% ����: ���������Ч����DDPʵ��
% ����: scale ---  ʵ����ģ��|C|��ֵ, 
%       max_gap --- CƬ����gap�Ͻ�
% ���������������ؼ���A,B,C,����DDP��Чʵ��

% �������������ʼ��
rand('twister',sum(100*clock));

% CƬ����gap�Ͻ�̶�Ϊmax_gap;
if max_gap<3
    max_gap=3;
end

SC=[0];SA=[0];SB=[0];
ja=2;jb=2;
if scale<2 
    scale=2 
end
LC=scale;
for jc=2:LC
    %disp('Next Gap:');
    gap=randint(1,max_gap);
    SC(jc)=SC(jc-1)+gap;
    r=rand;
    if r<0.45
        SA(ja)=SC(jc); ja=ja+1;
    elseif r<0.9
        SB(jb)=SC(jc);jb=jb+1;
    else
        SA(ja)=SC(jc);ja=ja+1;
        SB(jb)=SC(jc);jb=jb+1;
    end
end
%disp('Last Gap:');
gap=randint(1,max_gap);
SC(LC+1)=SC(LC)+gap;
SA(ja)=SC(LC+1);
SB(jb)=SC(LC+1);
A=[];B=[];C=[];
for i=1:length(SA)-1
    A(i)=SA(i+1)-SA(i);
end
for i=1:length(SB)-1
    B(i)=SB(i+1)-SB(i);
end
for i=1:length(SC)-1
    C(i)=SC(i+1)-SC(i);
end
%disp('New DDP Instance is generated: SA,SB,SC');
%SA
%SB
%SC
%disp('New DDP Instance is generated: unsorted A,B,C');
%A
%B
%C
% ����
A=sort(A);B=sort(B);C=sort(C);
%disp('New DDP Instance is generated: sorted A,B,C');
%A
%B
%C


% ���ʵ���ṹ
m=length(A);n=length(B);k=length(C);
error=0;
if m<2
    disp('Error: A is too short');
    error=1;
end
if n<2
    disp('Error: B is too short');
    error=1;
end
if k<2
    disp('Error: C is too short');
    error=1;
end

if error==1
    %��������
    disp('A, B or C is too short. Try again:');
    [A,B,C]=randDDPinstance(scale,max_gap);
end

% �޲����г���
m=length(A);n=length(B);k=length(C);
ZZ=[];
if k>m+n-1
    L=k-(m+n-1);
    ZZ(1:L)=0;
    A=[ZZ A];
    m=m+L;
else
    L=(m+n-1)-k;
    ZZ(1:L)=0;
    C=[ZZ C];
    k=k+L;
end

return
