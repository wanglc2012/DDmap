function [tt,g,f_opt,u,v,w]=permGA4(A,B,C,op)
% 
% ���ܣ����˫��������DDP(A,B,C)
% ���룺����������ؼ���A,B,C
%           �Ŵ��������ͼ��� op=[
%                          1 �û����Ͻ���
%                          2 �ο����򽻲�
%                          3 �������
%                          4 Ƭ�η�ת
%                          5 Ƭ����λ]
% ���������������ʱtt, ��������g, ������Ӧֵf_opt
%           ˫��������Ľ⣬�������û�(u,v,w)��ʾ
%

% ��ʼ��
tt=0;g=0;f_opt=0;u=[];v=[];w=[];

% ��������ʵ��
m=length(A);n=length(B);k=length(C);
if m<2
    disp('Error: A is too short');
    g=0;f_opt=0;u=[];v=[];
    return
end
if n<2
    disp('Error: B is too short');
    g=0;f_opt=0;u=[];v=[];
    return
end
if k<2
    disp('Error: C is too short');
    g=0;f_opt=0;u=[];v=[];
    return
end
% �޲����г���
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

% �趨�Ŵ����ʡ���Ⱥ��ģ������������
pc=0.85;
pm=2/(m+n);pm_gap=(0.55-pm)/200;
pop=50;maxG=10000;

% ��Ⱥ��ʼ��
% ������Ԫ������bu,bv��ʾ�䳤Ⱦɫ����Ⱥ
bu=cell(1,pop);
bv=cell(1,pop);
for i=1:pop
    %disp('Individual i:');    
    bu{i}=randperm(m);
    bv{i}=randperm(n);
end

% ������ʼ
tic
g=0;f_opt=0;i_opt=0;f(1:pop)=0;Loop=1;gf=[];
while (Loop==1)    
    for i=1:pop
        % Ⱦɫ�����
        u=bu{i};v=bv{i};
        % ��Ӧֵ����
        f(i)=DD_fitness(A,B,C,u,v);        
        % ����
        if f(i)>f_opt
            f_opt=f(i);
            i_opt=i;
        end
    end
    % ����    
    bu_opt=[];bu_opt=bu{i_opt};
    bv_opt=[];bv_opt=bv{i_opt};
    g=g+1;
    % ��ʾ��������ֵ
    %disp('--------- optimal fitness of current generation ---------');
    %str=sprintf('%d    %f', g, f_opt);
    %disp(str);
    gf(g)=f_opt;
    
    if (f_opt<1) && (g<maxG)       
        
        % �Ľ��Ķ���ѡ������,ֻ����תһ�μ������ȫ��ѡ��,������ȷ��ÿһ�������Ÿ��岻�ᱻ��̭
	    % ����ÿһ���屻ѡ��ĸ���
	    fp=f./sum(f);
	    % ������ּ��
        pop_choose=round(pop*0.8);
	    gap=1.0/pop_choose;
	    %������������ʼ��
	    toss=rand*gap;
	    j=0;jf=0.0;i=0;
	    while(i<pop_choose)
		    j=j+1;
		    if j>pop
			    % Խ��,��ת����
			    toss=rand*gap;
			    j=1;jf=0.0;
		    end
    	    jf=jf+fp(j);
    	    if toss<jf
			    % ѡ���j������
                i = i + 1;
                bu{i}=bu{j};bv{i}=bv{j};
        	    toss = toss + gap;
            end
        end
        % �����µ��������
        for i=pop_choose+1:pop-1
            bu{i}=randperm(m);
            bv{i}=randperm(n);
        end
        % ����
        bu{pop}=bu_opt;
        bv{pop}=bv_opt;
        
        % ִ�н���
        for i=1:2:pop_choose
            % Ⱦɫ�帱���Ʊ�
            xu=[];xv=[];yu=[];yv=[];
            xu=bu{i};xv=bv{i};
            yu=bu{i+1};yv=bv{i+1};            
           
            if ismember(1,op)                
            % �û����Ͻ���
            r=rand;
            if r<pc
                bu{i}=xu(yu);bu{i+1}=yu(xu);
                bv{i}=xv(yv);bv{i+1}=yv(xv);
            end
            end
            
            if ismember(2,op)
            % �ο����򽻲�
            r=rand;
            if r<pc
                j=randint(2,m-1);
                [xu,yu]=opPermCross(A,xu,yu,j);
                bu{i}=xu;bu{i+1}=yu;
            end
            r=rand;
            if r<pc
                j=randint(2,n-1);
                [xv,yv]=opPermCross(B,xv,yv,j);
                bv{i}=xv;bv{i+1}=yv;
            end
            end
        end
        
        % ִ�б���
        for i=1:pop_choose
            % Ⱦɫ�帱���Ʊ�
            xu=[];xv=[];
            xu=bu{i};xv=bv{i}; 
            
            if ismember(3,op)
            % �������
            r=rand;
            if r<pm
                j1=randint(1,m);j2=randint(1,m);
                t=xu(j1);xu(j1)=xu(j2);xu(j2)=t;
            end
            r=rand;
            if r<pm
                j1=randint(1,n);j2=randint(1,n);
                t=xv(j1);xv(j1)=xv(j2);xv(j2)=t;
            end  
            end            

            if ismember(4,op)
            % Ƭ�η�ת
            p=length(xu);q=length(xv);
            k2=randint(1,p);k1=randint(1,k2);L=k2-k1+1;
            r=rand;
            if r<pm
                xu(k1:1:k2)=xu(k2:-1:k1);
            end
            k2=randint(1,q);k1=randint(1,k2);L=k2-k1+1;
            r=rand;
            if r<pm
                xv(k1:1:k2)=xv(k2:-1:k1);
            end
            end
            
            if ismember(5,op)
            % Ƭ��ѭ����λ
            p=length(xu);q=length(xv);
            k2=randint(1,p);k1=randint(1,k2);L=k2-k1+1;
            r=rand;
            if r<pm
                t=xu(k1);
                xu(k1:k2-1)=xu(k1+1:k2);
                xu(k2)=t;
            end
            k2=randint(1,q);k1=randint(1,k2);L=k2-k1+1;
            r=rand;
            if r<pm
                t=xv(k1);
                xv(k1:k2-1)=xv(k1+1:k2);
                xv(k2)=t;
            end
            end

            % ���������
            bu{i}=xu;bv{i}=xv;            
        end
        % �����������
        pm=pm+pm_gap;
        if pm>=0.55
            pm=2/(m+n);
        end
    else
        Loop=0;
    end
end
tt=toc;
% ����Ⱦɫ�����
u=bu{i_opt};
v=bv{i_opt};
w=ddseq(A,B,C,u,v);
f_opt;
disp('[time gen fitness]=');
disp([tt, g, f_opt]);
return

function [f]=DD_fitness(A,B,C,u,v)
% ����(u,v)����Ӧֵ
% ��ʼ��
uA=[];vB=[];SuA=[];SvB=[];A_B=[];F=[];
% ���û�u�����ռ�A�е�Ԫ��
uA=A(u); 
% ���û�v�����ռ�B�е�Ԫ��
vB=B(v);
% �����ۼƼ���Sigma(u(A))
SuA=SigmaSet(uA);
% �����ۼƼ���Sigma(v(B))
SvB=SigmaSet(vB);
% ����ϲ�Sigma(u(A)),Sigma(v(B))
A_B=merging(SuA,SvB);
% ���㲽���
F=DeltaSet(A_B);
% ��F��������
F=sort(F);
% ����F��C��ƥ��̶�
[fitnum,missnum]=matching(F,C);
f=fitnum/(fitnum+missnum);
return

function [SA]=SigmaSet(A)
SA=[];
SA(1)=A(1);
for i=2:length(A)
    SA(i)=SA(i-1)+A(i);
end
return

function [DA]=DeltaSet(A)
DA=[];
DA(1)=A(1);
for i=2:length(A)
    DA(i)=A(i)-A(i-1);
end
return

function [T]=merging(A,B)
% �ϲ��������򼯺�A,B,��ȥ��ĩβ�ظ�Ԫ��
la=length(A);lb=length(B);
T=[];
i=1;j=1;k=1;
while (i<=la) || (j<=lb)
    if (i<=la) && (j<=lb)
        if A(i)<B(j)
            T(k)=A(i);i=i+1;k=k+1;
        else
            T(k)=B(j);j=j+1;k=k+1;
        end
    else
        while (i<=la)
            T(k)=A(i);i=i+1;k=k+1;
        end
        while (j<=lb)
            T(k)=B(j);j=j+1;k=k+1;
        end
    end
end
% ȥ��ĩβ�ظ�Ԫ��
T(length(T))=[];
return

function [fitnum,missnum]=matching(F,C)
% �����������򼯺�F,C��ƥ���ʧ���������
fitnum=0;missnum=0;
lf=length(F);lc=length(C);
i=1;j=1;
while (i<=lf) || (j<=lc)
    if (i<=lf) && (j<=lc)
        if F(i)<C(j)
            missnum=missnum+1;i=i+1;             
        elseif C(j)<F(i) 
            missnum=missnum+1;j=j+1;
        else
            fitnum=fitnum+1;i=i+1;j=j+1;        
        end
    else
        if (i<=lf)
            missnum=missnum+(lf-i+1);i=lf+1;
        elseif (j<=lc)
            missnum=missnum+(lc-j+1);j=lc+1;
        end
    end
end
return

function [w]=ddseq(A,B,C,u,v)
% ����˫�������м�, ����C��Ӧ���û�w
w=[];
uA=[];vB=[];wC=[];
%u
%v
uA=A(u);
vB=B(v);
SuA=[];SvB=[];
SuA=SigmaSet(uA);
SvB=SigmaSet(vB);
TT=[];
TT=merging(SuA,SvB);
DT=[];
DT=DeltaSet(TT);
[TT,iw]=sort(DT);
w=invperm(iw);
%C
wC=C(w);

% ��ӡ
disp('Double digest sequence:');
%disp('[u; A(u)]=');disp([u; uA]);
%disp('[v; B(v)]=');disp([v; vB]);
%disp('[w; C(w)]=');disp([w; wC]);
str=strABC(uA,vB,wC);
disp(str);
return

function [ix]=invperm(x)
ix=[];
lx=length(x);
ix(x)=1:lx;
return


