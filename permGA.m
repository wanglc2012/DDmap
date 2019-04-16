
function [tt,g,f_opt,u,v,w]=permGA(A,B,C,op)
% 
%Function: solve double digestion problem DDP(A,B,C)
%Input: three ordered multiple sets A,B, and C
%Set of genetic operator types op=[
%                          1 PCC
%                          2 RSC
%                          3 P4X
%                          4 FLP
%                          5 CSH]
% output: evolution computation time tt, evolution generations g, optimal fitness value f_opt
%The solution to the DDP ,it is represented by three permutation(u,v,w)
%




%  Initialization
tt=0;g=0;f_opt=0;u=[];v=[];w=[];

% op=[];
% op=input('please input op:' );
% op=[1,2,3,4,5,6,7];
% A =[ 1     1     2     2     2     4     4     5     5     6     6     7     8];
% B =[ 1     3     3     4     6    14    22];
% C =[1     1     1     1     1     2     2     2     2     2     2     3     3     4     4     5     5     6     6];



% Analyze input instances
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
% Repair sequence length
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

% Set genetic probability, population size and maximum evolution generation
pc=0.85;
pm=2/(m+n);pm_gap=(0.55-pm)/200;
pop=50;maxG=10000;

% InitializE Population 
% bu and bv were used to represent variable-length chromosome population
bu=cell(1,pop);
bv=cell(1,pop);
for i=1:pop
    %disp('Individual i:');    
    bu{i}=randperm(m);
    bv{i}=randperm(n);
end

% Evolution began
tic
g=0;f_opt=0;i_opt=0;f(1:pop)=0;Loop=1;gf=[];
while (Loop==1)    
    for i=1:pop
        % Chromosome decoding
        u=bu{i};v=bv{i};
        % Calculate the fitness value
        f(i)=DD_fitness(A,B,C,u,v);        
        % Select excellent individual
        if f(i)>f_opt
            f_opt=f(i);
            i_opt=i;
        end
    end
    % Preserve excellent individual
    bu_opt=[];bu_opt=bu{i_opt};
    bv_opt=[];bv_opt=bv{i_opt};
    g=g+1;
    % Shows contemporary optimal values
    %disp('--------- optimal fitness of current generation ---------');
    %str=sprintf('%d    %f', g, f_opt);
    %disp(str);
    gf(g)=f_opt;
   
  
    
    if (f_opt<1) && (g<maxG)       
        
        %Improved wheel selection operator, only need to rotate once to complete all the selection, and it ensures that each generation of the best individual will not be eliminated
        %Calculate the probability of each individual being selected
	    fp=f./sum(f);
	    % Calculate the wheel spacing
        pop_choose=round(pop*0.8);
	    gap=1.0/pop_choose;
	    %Calculate the random start point of the roulette wheel
	    toss=rand*gap;
	    j=0;jf=0.0;i=0;
	    while(i<pop_choose)
		    j=j+1;
		    if j>pop
			    %Over the border, turn the wheel again
			    toss=rand*gap;
			    j=1;jf=0.0;
		    end
    	    jf=jf+fp(j);
    	     if toss<jf
			    % Select the Jth individual
                i = i + 1;
                bu{i}=bu{j};bv{i}=bv{j};
        	    toss = toss + gap;
             end
        end
        % Introduce new random individuals
        for i=pop_choose+1:pop-1
            bu{i}=randperm(m);
            bv{i}=randperm(n);
        end
        % Preserve excellent individual
        bu{pop}=bu_opt;
        bv{pop}=bv_opt;
        
       
        % execute crossover
       
        
        
        for i=1:2:pop
            % Copy preparation
            xu=[];xv=[];yu=[];yv=[];
            xu=bu{i};xv=bv{i};
            yu=bu{i+1};yv=bv{i+1}; 
           
if ismember(1,op)                
            % PCC
%             a=1
            r=rand;
            if r<pc
                bu{i}=xu(yu);bu{i+1}=yu(xu);
                bv{i}=xv(yv);bv{i+1}=yv(xv);
            end
end

if ismember(2,op)
%      b=2
            % RSC
            r=rand;
            if r<pc
                j=randi([2,m-1]);
                [xu,yu]=opPermCross(A,xu,yu,j);
                bu{i}=xu;bu{i+1}=yu; 
            end
            r=rand;
            if r<pc
                j=randi([2,n-1]);
                [xv,yv]=opPermCross(B,xv,yv,j);
                bv{i}=xv;bv{i+1}=yv;
            end
end


if ismember(7,op)
 
L=length(xu);
Oi=zeros(1,L);
Oj=zeros(1,L);
Xi=xu;
Xj=yu;
for k=1:L
    r=rand;
    if Xi(k)==Xj(k)
        Oi(k)=Xi(k);
        Oj(k)=Xi(k);
    else if r<=0.5
         Oi(k)=Xj(k);
         Oj(k)=Xi(k);
            else 
            Oi(k)=Xi(k);
            Oj(k)=Xj(k);
        end
           end
end
end        
        
        
        end
            
        
        % execute muation
        for i=1:pop
            % Copy preparation
            xu=[];xv=[];
            xu=bu{i};xv=bv{i}; 
           
if ismember(3,op)
  
            % P4X
            r=rand;
            if r<pm
                j1=randi([1,m]);j2=randi([1,m]);
                t=xu(j1);xu(j1)=xu(j2);xu(j2)=t;
            end
            r=rand;
            if r<pm
                j1=randi([1,n]);j2=randi([1,n]);
                t=xv(j1);xv(j1)=xv(j2);xv(j2)=t;
            end  
end  

if ismember(4,op)
%     d=4
            % FLP
            p=length(xu);q=length(xv);
            k2=randi([1,p]);k1=randi([1,k2]);L=k2-k1+1;
            r=rand;
            if r<pm
                xu(k1:1:k2)=xu(k2:-1:k1);
            end
            k2=randi([1,q]);k1=randi([1,k2]);L=k2-k1+1;
            r=rand;
            if r<pm
                xv(k1:1:k2)=xv(k2:-1:k1);
            end
            
end



if ismember(5,op)
%     e=5
            % CSH
            p=length(xu);q=length(xv);
            k2=randi([1,p]);k1=randi([1,k2]);L=k2-k1+1;
            r=rand;
            if r<pm
                t=xu(k1);
                xu(k1:k2-1)=xu(k1+1:k2);
                xu(k2)=t;
                ma=1;
            end
           k2=randi([1,q]);k1=randi([1,k2]);L=k2-k1+1;
            r=rand;
            if r>pm
                t=xv(k1);
                xv(k1:k2-1)=xv(k1+1:k2);
                xv(k2)=t;
            end
       
end

 

if ismember(6,op)
%     f=6
      % 2005 S-K muation
      A=xu;
     for i=1:m 
            r=rand;
            if r<pm
                j=randi([1,m]);
                t=A(i); A(i)=A(j);A(j)=t;
            end
            pm_gap=(0.45-pm)/200;
            pm=pm+pm_gap;
        if pm>=0.45
            pm=pm-pm_gap;
            if pm<2/(m+n)
               pm=2/(m+n); 
            end
        end
     end
end


            % Preserve muation results
            bu{i}=xu;bv{i}=xv;            
        end
        % Modified muation probability 
        pm=pm+pm_gap;
        if pm>=0.55
            pm=2/(m+n);
        end
    else
        Loop=0;

end
end
    

tt=toc;
% Optimal chromosome decoding
u=bu{i_opt};
v=bv{i_opt};
w=ddseq(A,B,C,u,v);
f_opt;
disp('[time gen fitness]=');
disp([tt, g, f_opt]);
return



function [f]=DD_fitness(A,B,C,u,v)
% Calculate the fitness value of (u,v)
% initialization
uA=[];vB=[];SuA=[];SvB=[];A_B=[];F=[];
% Recollect elements of A by permutation u
uA=A(u); 
% Recollect elements of B by permutation v
vB=B(v);
% Compute cumulative set Sigma(u(A))
SuA=SigmaSet(uA);
% Compute cumulative setSigma(v(B))
SvB=SigmaSet(vB);
% Sigma(u(A)),Sigma(v(B))
A_B=merging(SuA,SvB);
% Calculate the step difference set
F=DeltaSet(A_B);
% Sort F
F=sort(F);
% Calculate the match between F and C
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
% Merge two ordered sets A,B, and remove the end repeating elements
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
% Remove the trailing repeating element
T(length(T))=[];
return

function [fitnum,missnum]=matching(F,C)
% The number of matching and mismatch components of two ordered sets F and C is calculated
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
% Calculate the double digestion sequence and return the replacement w corresponding to C
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

% ´òÓ¡
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


