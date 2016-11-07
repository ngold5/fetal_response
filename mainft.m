clear all
close all
global Hbm c1m c2m Vut Bput Rut
global T0 fv0 fvh tuco tucor ttr
% initial values


Y0=[52.3;11.5;89.4;11.5;11.5;...
    14.6;17.6;11.50;12.;0.33;...
    29.5;0.043;45.5;6.3;-10.5;...
    10.;52.5;111;7.5;196;...
    24.6;0.54;0.57;0.61;0.63;...
    13;1;2.1;0.12;47];

%{
Y0=[60;14;70;12;14;...
    16.0;16;18.;12;.4;...
    30;0.01;40;5;-5.;...
    10.5;60;90;10;40;...
    15;0.52;0.52;0.52;0.52;...
    8; 1.2;1.8;0.12;50;...
    1.2;1.8;0.12;50];
%}

%Y0=[60;14;70;12;14;...
%    16.0;16;.06;15.;.4;...
%    30.];
% add oxygen distribution, now the variables are partial pressure
alf = 1.34;
bet = 0.031*0.1;
Hbf = 1.7e-1;
Hbm = 1.2e-1;
c1f = 1.04e4;
c2f = 150;
c1m = 2.34e4;
c2m = 150;
Vut = 500;
Bput = 76;
Rut = 0.07*60*1;

% regulation para
Dtv = 0.2;  % s  delay time
Gtv = 0.09;  % s^2
Gtv = 0.04;
tauv = 1.5;  % s
taus = 3;
Dts = -0.2;  % s
Gts = -0.13;

T0 = 0.21;  %s  T = T0+DTs+DTv
%T0 = 0.16;
%T0 = 0.4;
fv0 = 5.85;  %1/s
%fv0 = 2.50;
fvh = 5;  % 1/s

tuco = 60*1;
tucor = 60*1.5;  % recovery time
%tuco = 100;
%tucor = 50;  % recovery time
ttr = 600;
%ttr = 0;
%T = ttr + 60*20.;
T = ttr + (tuco+tucor)*30;
T = ttr + 60*60*4.5;
%T = 600;
dt = 2.0;
tspan = 1:dt:T;

% parameters
factor = 1/133.322;
factor2 = 1000/133.322;

u(1) = 6*factor;  % Ra
u(2) = 30*factor;  % Rv
u(3) = 200*factor;  % Rmc0
u(4) = 698*factor;   % Rummc
u(5) = 2*factor;  % Rumv
u(6) = 1398*factor;  % Rcmc0
u(7) = 2*factor;  % Rcv
u(8) = 2/factor2*4;   % Ca
u(9) = 35/factor2;    % Cv
u(10) = 11/factor2;   % Cum
u(11) = 0.57/factor2; % Cc
% regulation
u(12) = 3.2;  % fevmin
u(13) = 6.3;  % fevmax
u(14) = 0.2;   % Wcv
u(15) = 1;      % Wbv
u(16) = 7.06;    % kev
u(17) = 1.7;   % kvh
u(18) = 4.09;   % kab
u(19) = 7.24;   % kac
u(20) = 25;   %fabn
u(21) = 11;   %facn
u(22) = 0.0675;  % kes
u(23) = 1;  % wby
u(24) = 1;   % Wcy
u(25) = 60;   % fesmax
u(26) = 2.1;   % fesinf
u(27) = 16.11;   % fes0
u(28) = 2.66;   % fesmin
u(29) = 7;   % pO20
u(30) = 0.018;   % kcmc
u(31) = 0.05;    % cOan
u(32) = 9;  % fspn
u(33) = 1.44;   % kkmc
u(34) = 40;   % pn
u(35) = 7.5;   % pO2n
u(36) = 10;    % fvn
u(37) = 5.5;   % fshn
u(38) = Gts;
u(39) = Gtv;
% metabolic
u(40) = c1f;
u(41) = c2f;
u(42) = alf;
u(43) = Hbf;
u(44) = bet;

u(45) = 9.33;  % Kf  for O2 met
u(46) = 0.068;  % cOth
u(47) = 0.31;   % Omet0

u(48) = 0.25;   % klc
u(49) = 50;   % kk5
u(50) = 0.04;  % kk1
u(51) = 0.02;  % kk2
u(52) = 0.01;  % kk4
u(53) = 30;  % kk7

u(54) = 0.042;   % Omet0c
u(55) = 0.039;    % cOthc
u(56) = 1;  % KK coef
% supplement

u(57) = tauv;
u(58) = taus;
u(59) = 30;  % tauisc
u(60) = 0.05;  % psi
%
u(61) = 1;  % glc0
u(62) = .8;   %lac0
u(63) = 0.12;  % py0
u(64) = 40;  % hc0

%
u(65) = 10;  % pyruvate consumption
u(66) = 44;  % CO2mc0
u(67) = 0.1;  % H+ feedback to CO2
u(68) = 0.1;

u = u';

tic
options = odeset('MaxStep',1e3,'RelTol',1e-4);
[t,Y] = ode15s(@(t,Y)mainode1(t,Y,u),tspan,Y0,options);
%[t,Y,dYdU] = sens_ind('mainode_sens',tspan,Y0,[],u);
%lags = [0.2, 2]; 
%lags = [];

%tic 

%sol = dde23(@maineqns81dde,lags,Y0,tspan);

toc
% get O2 concentration in umbilical cord

%t = sol.x;
%Y = sol.y;
%Y = Y';

figure(1)
subplot(4,1,1)
plot(t,Y(:,1),'r')
%xlim([ttr T])
subplot(4,1,2)
plot(t,Y(:,10),'b')
%xlim([ttr T])
subplot(4,1,3)
plot(t,Y(:,12),'b')
%xlim([ttr T])
subplot(4,1,4)
plot(t,Y(:,9),'b')
%xlim([ttr T])
%xlim([10 20])
print -depsc ABP

figure(2)
subplot(4,1,1)
plot(t,Y(:,6),'b')
%xlim([ttr T])
subplot(4,1,2)
plot(t,Y(:,7))
subplot(4,1,3)
plot(t,Y(:,8),'b')
%xlim([ttr T])
subplot(4,1,4)
plot(t,Y(:,9),'b')
%xlim([ttr T])
print -deps MpO2


kab = 4.09;
fabmin = 2.52;
fabmax = 47.78;
fabn = 40;
Gab = (Y(:,13)-fabn)./(Y(:,17)-10);
%kab = (fabmax-fabmin)/4./Gab;
fabs = (fabmin + fabmax*exp((Y(:,17)-15)./kab))./(1+exp((Y(:,17)-15)./kab));

figure(3)
%FHR = 1./(T0-0.2+Y(:,10));
%FHR = 1./(T0-0.177+Y(:,12)+Y(:,10));
FHR = 1./(T0-0.13*log(5.5-2.66+1)+Y(:,12)+Y(:,10));
subplot(2,1,1)
plot(t,FHR*60) % heart rate  bpm
%xlim([ttr T])
ylabel FHR
subplot(2,1,2)
plot(t,Y(:,17))
%plot(t,fabs)
print -deps ffhr

m = round(ttr/dt);
L = length(t);
time = t(1:L);
for i = 1 : length(time)
sw(i)=1;
if time(i)>ttr && ...
   (time(i)-floor(time(i)/(tuco+tucor))*(tuco+tucor))<=tuco
    sw(i) = 0;
end
end

%gefun = Ge*0.985+0.005*Ge*sw;
%plot(time-20,gefun)
Sa = 1 - c1f./(c1f + Y(:,6).^3+c2f*Y(:,6));
cOa = alf*Hbf*Sa + bet*Y(:,6);  % 

Sc = 1 - c1f./(c1f + Y(:,9).^3+c2f*Y(:,9));
cOc = alf*Hbf*Sc + bet*Y(:,9);  % c
Sum = 1 - c1f./(c1f + Y(:,7).^3+c2f*Y(:,7));
cOum = alf*Hbf*Sum + bet*Y(:,7);  % um
figure(4)
subplot(4,1,1)
plot(t,cOc,'--')  % cerebral O2
hold on
plot(t,cOum,'-')
plot(t,cOa,'r-')
plot(t,Y(:,8),'m-')
%xlim([ttr-50 T])
subplot(4,1,2)
plot(t,Y(:,1))
%xlim([ttr-50,T])
subplot(4,1,3)
plot(t,FHR*60)
%xlim([ttr-50 T])
subplot(4,1,4)
%plot(time+tuco,sw)
plot(time,sw)
%axis([ttr-50 T -0.1 1.1])

print -depsc datashow

figure(5)
subplot(4,1,1)
plot(t,Y(:,12),'r')
%xlim([ttr T])
subplot(4,1,2)
plot(t,Y(:,13),'b')
%xlim([ttr T])
subplot(4,1,3)
plot(t,Y(:,14),'b')
%xlim([ttr T])
subplot(4,1,4)
plot(t,Y(:,16),'b')
%xlim([ttr T])%xlim([10 20])
print -depsc regulation

figure(6)
subplot(4,1,1)
plot(t,Y(:,18),'r')
subplot(4,1,2)
plot(t,Y(:,19),'r')
subplot(4,1,3)
plot(t,Y(:,20),'r')
subplot(4,1,4)
plot(t,Y(:,21),'r')
print -depsc motherP

figure(7)
subplot(4,1,1)
plot(t,Y(:,22))
%xlim([ttr T])
subplot(4,1,2)
plot(t,Y(:,23))
%xlim([ttr T])
subplot(4,1,3)
plot(t,Y(:,24))
%xlim([ttr T])
subplot(4,1,4)
plot(t,Y(:,25))
%xlim([ttr T])
print -depsc metab1

figure(8)
subplot(4,1,1)
plot(t,Y(:,27))
%xlim([ttr T])
subplot(4,1,2)
plot(t,Y(:,28))
%xlim([ttr T])
subplot(4,1,3)
plot(t,Y(:,29))
%xlim([ttr T])
subplot(4,1,4)
plot(t,Y(:,30))
%xlim([ttr T])
print -depsc metab2

% calculate BE
pco2 = (Y(:,24) - 0.244)/0.007;
pH = -log(Y(:,30)*1e-9);
BE = 0.02786*pco2.*10.^(pH-6.1)+13.77*pH - 124.58;
figure(10)
plot(t,BE)
print -depsc BDval

fid = fopen('time.txt', 'w');
fprintf(fid,'% 6.14f ', t(:));
fclose(fid);

fid = fopen('arterialpr.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,1));
fclose(fid);
%fid = fopen('Pv.txt', 'w');
%fprintf(fid,'% 6.14f ', Y(:,2));
%fclose(fid);
%fid = fopen('Vv.txt', 'w');
%fprintf(fid,'% 6.14f ', Y(:,3));
%fclose(fid);
%fid = fopen('umpr.txt', 'w');
%fprintf(fid,'% 6.14f ', Y(:,4));
%fclose(fid);
%fid = fopen('cpr.txt', 'w');
%fprintf(fid,'% 6.14f ', Y(:,5));
%fclose(fid);
fid = fopen('pcOa.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,6));
fclose(fid);
fid = fopen('pcOum.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,7));
fclose(fid);
fid = fopen('pcOmc.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,8));
fclose(fid);
fid = fopen('pcOc.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,9));
fclose(fid);
fid = fopen('DTv.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,10));
fclose(fid);
%fid = fopen('PO2m.txt', 'w');
%fprintf(fid,'% 6.14f ', Y(:,11));
%fclose(fid);
fid = fopen('DTs.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,12));
fclose(fid);
%fid = fopen('fabs.txt', 'w');
%fprintf(fid,'% 6.14f ', Y(:,13));
%fclose(fid);
%fid = fopen('facs.txt', 'w');
%fprintf(fid,'% 6.14f ', Y(:,14));
%fclose(fid);
%fid = fopen('fsh0.txt', 'w');
%fprintf(fid,'% 6.14f ', Y(:,15));
%fclose(fid);
%fid = fopen('Rcmc.txt', 'w');
%fprintf(fid,'% 6.14f ', Y(:,16));
%fclose(fid);
fid = fopen('Pabar.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,17));
fclose(fid);
%fid = fopen('Pam.txt', 'w');
%fprintf(fid,'% 6.14f ', Y(:,18));
%fclose(fid);
%fid = fopen('Pvm.txt', 'w');
%fprintf(fid,'% 6.14f ', Y(:,19));
%fclose(fid);
%fid = fopen('Vlm.txt', 'w');
%fprintf(fid,'% 6.14f ', Y(:,20));
%fclose(fid);
%fid = fopen('Pivs.txt', 'w');
%fprintf(fid,'% 6.14f ', Y(:,21));
%fclose(fid);
fid = fopen('CO2a.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,22));
fclose(fid);
fid = fopen('CO2um.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,23));
fclose(fid);
fid = fopen('CO2mc.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,24));
fclose(fid);
fid = fopen('CO2c.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,25));
fclose(fid);
fid = fopen('glucose.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,27));
fclose(fid);
fid = fopen('lac.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,28));
fclose(fid);
fid = fopen('pyr.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,29));
fclose(fid);
fid = fopen('Hydc.txt', 'w');
fprintf(fid,'% 6.14f ', Y(:,30));
fclose(fid);

fid = fopen('ucoinfo.txt', 'w');
fprintf(fid,'% 6.14f ', sw(:));
fclose(fid);

% sensitivity
% plot(t,dYdU(:,j,k)./Y(:,j)*u(k))


%fid = fopen('Smap2.txt', 'w');
%fprintf(fid,'% 6.14f ', dYdU(:,17,2));
%fclose(fid);
% plot(t,dYdU(:,j,k)./Y(:,j)*u(k))
senst = [];
%for jj=1:length(u)
%    senst(jj,1:30) = max(dYdU(:,:,jj));
%end
for jj=1:length(u)
    for ii = 1:30
        tmp(ii,jj) = max(dYdU(:,ii,jj)./Y(:,ii)*u(jj));
    end
    senst(jj) = max(tmp(:,jj));
%    senst(jj) = max(dYdU(:,:,jj));
end
save('senranking0.txt','senst','-ASCII')
fid = fopen('Senranking.txt', 'wt');
for i=1:length(u)
  fprintf(fid,'% 6.14f ', senst(i,:));
end
fclose(fid);





