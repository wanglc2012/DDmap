function [str]=strABC(A,B,C)
la=length(A);
lb=length(B);
lc=length(C);
str=[];
str=strcat(str, 'A:=[',  sprintf('%d, ', A(1:la-1)), sprintf('%d', A(la)),']: ');
str=strcat(str, 'B:=[',  sprintf('%d, ', B(1:lb-1)), sprintf('%d', B(lb)),']: ');
str=strcat(str, 'C:=[',  sprintf('%d, ', C(1:lc-1)), sprintf('%d', C(lc)),']: ');
%disp(str);
return 



