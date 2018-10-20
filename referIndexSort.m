function [IA,IAI]=referIndexSort(IA,A,B)
% 参照索引排序: 参照B,对A的索引IA（未必是连续的）引进行排序，要求结果与B中的相容，即如果有A(IA(i))<A(IA(j))，则一定能B中找到这样的顺序关系
%              必要条件length(IA)<=length(A)
% 输出：排序后的索引IA和“索引IA的索引IAI”

L=length(IA);IAI=1:L;

% 参考索引冒泡排序
swap=1;
while (swap==1)
    swap=0;
for i=1:L-1
    for j=i+1:L
        if refOrder(A(IA(IAI(j-1))),A(IA(IAI(j))),B)>0
            t=IAI(j-1);IAI(j-1)=IAI(j);IAI(j)=t;
            swap=1;
        end
    end
end
end

% % 改进的冒泡参考索引排序
% n = LIA;IAI=1:n;
% while (n>0)
%     newn = 0;
%     for j = 2:n
%         if refOrder(A(IA(IAI(j-1))),A(IA(IAI(j))),B)>0
%             t=IAI(j-1);IAI(j-1)=IAI(j);IAI(j)=t;
%             newn = j;
%         end
%     end
%     n = newn;
% end

IA=IA(IAI);
return

function [k]=refOrder(x,y,B)
% if x<=y
%     k=0;
% else
%     k=1;
% end
k=0;LB=length(B);
jx=1;
while (x~=B(jx) && jx<LB) 
    jx=jx+1;
end
if x~=B(jx)    
    return
end
jy=LB;
while (y~=B(jy) && jy>1) 
    jy=jy-1;
end
if y~=B(jy)
    return
end
if jx>jy
    k=1;
end
return

       