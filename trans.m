%Scaling-rounding-adjusting approach for Cases of INS2

a = [5509 5626 6527 6766 7233 16841];
b = [3526 4878 5643 5804 7421 21230];
c = [1120 1868 2564 2752 3240 3526 3758 3775 4669 5509 15721];
aa = a*0.001;
bb = b*0.001;
cc = c*0.001;
for i=0:5
    k = 0.1*i;
    ar = roundn(aa+k,0);
    br = roundn(bb+k,0);
    cr = roundn(cc+k,0);
    as = sum(ar);
    bs = sum(br);
    cs = sum(cr);
    if as==bs && bs==cs
        ar,br,cr
        disp('Get it!')
        break;
    end;
   
end; 
for i=0:5
    k = 0.1*i;
    ar = roundn(aa - k,0);
    br = roundn(bb - k,0);
    cr = roundn(cc - k,0);
    as = sum(ar);
    bs = sum(br);
    cs = sum(cr);
    if as==bs && bs==cs
        ar,br,cr
        disp('Get it!')
        break;
    end;
    
end;    

A=ar;B=br;C=cr;
op=input('Please input the operator you want, op£º');
[tt,g,f_opt,u,v,w]=permGA(A,B,C,op);
A=a(u);B=b(v);C=c(w);
A
B
C




