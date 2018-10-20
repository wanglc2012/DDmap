function [newx,newy]=opPermCross(M,x,y,j)
% 保序置换交叉
% 输入：M --- 染色体基因座定义母版
%       x,y --- 两个置换，用来定义参与交叉的两条染色体 M(x) 和 M(y)
%       j ---  交叉点
% 输出：newx,newy --- 两个新的置换，由 M(newx),M(newy)给出，满足如下交叉效果
%                     M(newx)在j点之前继承自M(x),在j点之后参照M(y)进行了排序
%                     M(newy)在j点之前继承自M(y),在j点之后参照M(x)进行了排序

L=length(M);A=M(x);B=M(y);
newx=[];newy=[];
newx=referIndexSort(x(j+1:L),M,B);
newy=referIndexSort(y(j+1:L),M,A);
newx=[x(1:j) newx];
newy=[y(1:j) newy];
return








