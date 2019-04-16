%Comparison of DDmap and other algorithms under the condition of INS2

tt =[13.606   0.561   7.297   7.265   10.835   191.851   43.562];

gg =[7521.10   236.70   3934.30   3543.70   5315.60   100000   23967.20];

ss =[ 1.0   1.0  1.0  1.0   1.0  0   0.90];

x=[1 2 3 4 5 6 7 ];
figure;
% subplot(1,2,1);
plot(x,tt,'*-') ;
xlabel('Genetic operators');
ylabel('Average Running Time/(s) ');
% legend('INS2');
% title('(a)');

figure;
% subplot(1,2,2);
plot(x,ss,'*-') ;
xlabel('Genetic operators');
ylabel('Success Rate ');
% legend('INS2');
% title('(b)');
