function [tt,g,f_opt,u,v,w]=permGA4(A,B,C,op)
% 
% 功能：求解双消化问题DDP(A,B,C)
% 输入：三个有序多重集合A,B,C
%           遗传算子类型集合 op=[
%                          1 置换复合交叉
%                          2 参考保序交叉
%                          3 单点变异
%                          4 片段翻转
%                          5 片段移位]
% 输出：进化计算用时tt, 进化代数g, 最优适应值f_opt
%           双消化问题的解，用三个置换(u,v,w)表示
%

% 初始化
tt=0;g=0;f_opt=0;u=[];v=[];w=[];

% 解析输入实例
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
% 修补序列长度
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

% 设定遗传概率、种群规模和最大进化代数
pc=0.85;
pm=2/(m+n);pm_gap=(0.55-pm)/200;
pop=50;maxG=10000;

% 种群初始化
% 用两个元胞数组bu,bv表示变长染色体种群
bu=cell(1,pop);
bv=cell(1,pop);
for i=1:pop
    %disp('Individual i:');    
    bu{i}=randperm(m);
    bv{i}=randperm(n);
end

% 进化开始
tic
g=0;f_opt=0;i_opt=0;f(1:pop)=0;Loop=1;gf=[];
while (Loop==1)    
    for i=1:pop
        % 染色体解码
        u=bu{i};v=bv{i};
        % 适应值计算
        f(i)=DD_fitness(A,B,C,u,v);        
        % 择优
        if f(i)>f_opt
            f_opt=f(i);
            i_opt=i;
        end
    end
    % 保优    
    bu_opt=[];bu_opt=bu{i_opt};
    bv_opt=[];bv_opt=bv{i_opt};
    g=g+1;
    % 显示当代最优值
    %disp('--------- optimal fitness of current generation ---------');
    %str=sprintf('%d    %f', g, f_opt);
    %disp(str);
    gf(g)=f_opt;
    
    if (f_opt<1) && (g<maxG)       
        
        % 改进的赌轮选择算子,只需旋转一次即可完成全部选择,而且它确保每一代的最优个体不会被淘汰
	    % 计算每一个体被选择的概率
	    fp=f./sum(f);
	    % 计算赌轮间距
        pop_choose=round(pop*0.8);
	    gap=1.0/pop_choose;
	    %计算赌轮随机起始点
	    toss=rand*gap;
	    j=0;jf=0.0;i=0;
	    while(i<pop_choose)
		    j=j+1;
		    if j>pop
			    % 越界,重转赌轮
			    toss=rand*gap;
			    j=1;jf=0.0;
		    end
    	    jf=jf+fp(j);
    	    if toss<jf
			    % 选择第j个个体
                i = i + 1;
                bu{i}=bu{j};bv{i}=bv{j};
        	    toss = toss + gap;
            end
        end
        % 引入新的随机个体
        for i=pop_choose+1:pop-1
            bu{i}=randperm(m);
            bv{i}=randperm(n);
        end
        % 保优
        bu{pop}=bu_opt;
        bv{pop}=bv_opt;
        
        % 执行交叉
        for i=1:2:pop_choose
            % 染色体副本制备
            xu=[];xv=[];yu=[];yv=[];
            xu=bu{i};xv=bv{i};
            yu=bu{i+1};yv=bv{i+1};            
           
            if ismember(1,op)                
            % 置换复合交叉
            r=rand;
            if r<pc
                bu{i}=xu(yu);bu{i+1}=yu(xu);
                bv{i}=xv(yv);bv{i+1}=yv(xv);
            end
            end
            
            if ismember(2,op)
            % 参考保序交叉
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
        
        % 执行变异
        for i=1:pop_choose
            % 染色体副本制备
            xu=[];xv=[];
            xu=bu{i};xv=bv{i}; 
            
            if ismember(3,op)
            % 单点变异
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
            % 片段翻转
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
            % 片段循环移位
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

            % 保存变异结果
            bu{i}=xu;bv{i}=xv;            
        end
        % 修正变异概率
        pm=pm+pm_gap;
        if pm>=0.55
            pm=2/(m+n);
        end
    else
        Loop=0;
    end
end
tt=toc;
% 最优染色体解码
u=bu{i_opt};
v=bv{i_opt};
w=ddseq(A,B,C,u,v);
f_opt;
disp('[time gen fitness]=');
disp([tt, g, f_opt]);
return

function [f]=DD_fitness(A,B,C,u,v)
% 计算(u,v)的适应值
% 初始化
uA=[];vB=[];SuA=[];SvB=[];A_B=[];F=[];
% 按置换u重新收集A中的元素
uA=A(u); 
% 按置换v重新收集B中的元素
vB=B(v);
% 计算累计集合Sigma(u(A))
SuA=SigmaSet(uA);
% 计算累计集合Sigma(v(B))
SvB=SigmaSet(vB);
% 按序合并Sigma(u(A)),Sigma(v(B))
A_B=merging(SuA,SvB);
% 计算步进差集
F=DeltaSet(A_B);
% 对F进行排序
F=sort(F);
% 计算F和C的匹配程度
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
% 合并两个有序集合A,B,并去掉末尾重复元素
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
% 去掉末尾重复元素
T(length(T))=[];
return

function [fitnum,missnum]=matching(F,C)
% 计算两个有序集合F,C的匹配和失配分量个数
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
% 计算双消化序列及, 返回C对应的置换w
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

% 打印
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


