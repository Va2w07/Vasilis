refd = ref3d;
reft =  ref3t;
cold = col3d;
colt = col3t;
cont = foc3t;
cond = foc3d;
cont2 = foc3t2;
cond2 = foc3d2;

refd(isnan(refd)) = 0;
cold(isnan(cold)) = 0;
cond(isnan(cond)) = 0;
cond2(isnan(cond2)) = 0;

reft = reft - reft(1);
colt = colt - reft(1);
cont = cont - reft(1);
cont2 = cont2 - reft(1);

dt = (reft(351)-reft(1))./351;

reft = (reft(1):dt:(dt*351*2-1+reft(1)))';
colt = (colt(1):dt:(dt*351*2-1+colt(1)))';
cont = (cont(1):dt:(dt*351*2-1+cont(1)))';
cont2 = (cont2(1):dt:(dt*351*2-1+cont2(1)))';

refd = refd - refd(1);
cold = cold - cold(1);
cond = cond - cond(1);
cond2 = cond2 - cond2(1);

refd = padarray(refd,351,0,'post');
cold = padarray(cold,351,0,'post');
cond = padarray(cond,351,0,'post');
cond2 = padarray(cond2,351,0,'post');

csvwrite('ref3.dat', [reft,refd]);
csvwrite('col3.dat', [colt,cold]);
csvwrite('foc3.dat', [cont,cond]);
csvwrite('foc32.dat', [cont2,cond2]);
