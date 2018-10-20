
clear
clc
tt=[];gg=[];ss=[];
for i=10:10:100
[t,g,s]=simu1007(10,i);
tt=[tt, t];
gg=[gg, g];
ss=[ss,s];
end

figure;
X=10:10:100;
plot(X,tt,'b-o');
xlabel('Maximum Length of Piece in Sequence C');
ylabel('Average Running Time of Finding Exact DDP Solutions (sec)');
grid;

figure;
plot(X,gg,'r-d');
xlabel('Maximum Length of Piece in Sequence C');
ylabel('Average Evoluting Generations of Finding Exact DDP Solutions');
grid;

figure;
plot(X,ss,'r-d');
xlabel('Maximum Length of Piece in Sequence C');
ylabel('Success Rate of Finding Exact DDP Solutions');
grid;



