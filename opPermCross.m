function [newx,newy]=opPermCross(M,x,y,j)
% �����û�����
% ���룺M --- Ⱦɫ�����������ĸ��
%       x,y --- �����û�������������뽻�������Ⱦɫ�� M(x) �� M(y)
%       j ---  �����
% �����newx,newy --- �����µ��û����� M(newx),M(newy)�������������½���Ч��
%                     M(newx)��j��֮ǰ�̳���M(x),��j��֮�����M(y)����������
%                     M(newy)��j��֮ǰ�̳���M(y),��j��֮�����M(x)����������

L=length(M);A=M(x);B=M(y);
newx=[];newy=[];
newx=referIndexSort(x(j+1:L),M,B);
newy=referIndexSort(y(j+1:L),M,A);
newx=[x(1:j) newx];
newy=[y(1:j) newy];
return








