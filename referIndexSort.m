function [IA,IAI]=referIndexSort(IA,A,B)
% ������������: ����B,��A������IA��δ���������ģ�����������Ҫ������B�е����ݣ��������A(IA(i))<A(IA(j))����һ����B���ҵ�������˳���ϵ
%              ��Ҫ����length(IA)<=length(A)
% ���������������IA�͡�����IA������IAI��

L=length(IA);IAI=1:L;

% �ο�����ð������
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

% % �Ľ���ð�ݲο���������
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

       