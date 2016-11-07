clear all
close all
global Hbm c1m c2m Vut Bput Rut
global T0 fv0 fvh tuco tucor ttr

% set up initial values

Y0 = [52.3; 11.5; 89.4; 11.5; 11.5;...
	14.6; 17.6; 11.5; 12.0; 0.33;...
	29.5; 0.043; 45.5; 6.3; -10.5;...
	10; 52.5; 111; 7.5; 196;...
	24.6; 0.54; 0.57; 0.61; 0.62;...
	13; 1; 2.1; 0.12; 47];

% add in oxygen distribution
% the variables are partial pressures - volume included in units
alf = 1.34;
bet = 0.031*0.1;
Hbf = 1.7e-1;		% Couto
Hbm = 1.2e-1;		% Couto
c1f = 1.04e4;		% Couto
c2f = 150;			% Couto
c1m = 2.34e4;		% Couto
c2m = 150;			% Couto
Vut = 500;			
Bput = 76;
Rut = 0.07 * 60;

% regulation parameters
Dtv = 0.2;		% s 	Delta Tv - vagal contribution to heart period
Gtv = 0.09;		% s^2	gain
Gtv = 0.04;
tauv = 1.4;		% s
taus = 3;		% s
Dts = -0.2;		% s 	Delta Ts - sympathetic contribution to heart period
Gts = -0.13;	% s^2	gain

T0 = 0.21;		% s 	T = T0 + DTs + DTv
fv0 = 5.85;		% 1/s
fvh = 5;		% 1/s
	
tuco = 60;			% UCO length in seconds
tucor = 60*1.5;		% recovery time
ttr = 600;

% set up timing information now for simulation
T = ttr + (tuco + tucor) * 30;
T = ttr + 3600*4.5;

% time step
dt = 2.0;
tspan = 1:dt:T;

% ----------------% ----------------% ----------------% ----------------
% parameters
factor = 1 / 133.322;
factor2 = 1000 / 133.322;

% define the holder u for all of the parameters and conditions now
% this is called uu in the function handle call of qimingode.m
% the script that contains the ODE

% circulation component
u(1) = 6*factor;		% Ra 	arterial resistance
u(2) = 30*factor;		% Rv 	venous resistance
u(3) = 200*factor;		% Rmc0  systemic resistance
u(4) = 698*factor;		% Rummc umbilical cord resistance
u(5) = 2*factor;		% Rumv  umbilical venous resistance
u(6) = 1398*factor;		% Rcmc0 cerebral resistance
u(7) = 2*factor;		% Rcv 	cerebral venous resitance
u(8) = 2/factor2 * 4;	% Ca 	arterial compliance
u(9) = 35/factor2;		% Cv	venous compliance
u(10) = 11/factor2;		% Cum	umbilical compliance
u(11) = 0.57/factor2;	% Cc 	cerebral compliance

% regulation component
u(12) = 3.2;		% fevmin	vagal firing rate minimum
u(13) = 6.3;		% fevmax	vagal firing rate maximum
u(14) = 0.2;		% Wcv 		chemoreceptor weight parameter
u(15) = 1;			% Wbv 		baroreceptor weight parameter
u(16) = 7.06;		% kev		slope parameter in sigmoidal curve
u(17) = 1.7;		% kvh 		vagal slope parameter
u(18) = 4.09;		% kab 		baroreceptor slope for sigmoid
u(19) = 7.24; 		% kac		chemoreceptor slope for sigmoid
u(20) = 25;			% fabn		baroreceptor reference firing rate
u(21) = 11;			% facn		chemoreceptor reference firing rate
u(22) = 0.0675;		% kes		slope paramter for sigmoid
u(23) = 1;			% wby 		weight - baro
u(24) = 1;			% wcy		weight - chemo
u(25) = 60;			% fesmax	sympathetic firing rate maximum
u(26) = 2.1;		% fesinf	sympathetic firing rate parameter
u(27) = 16.11;		% fes0		sympathetic firing rate reference
u(28) = 2.66;		% fesmin	sympathetic firing rate minimum
u(29) = 7;			% pO20		initial oxygen partial pressure
U(30) = 0.018;		% kcmc		autoregulation resistance slope in sigmoid
u(31) = 0.05;		% cOan 		arterial oxygen concentration reference
u(32) = 9;			% fspn		alpha-sympth reference firing rate
u(33) = 1.44;		% kkmc		systemic slope in sigmoid
u(34) = 40;			% pn		reference partial pressure
u(35) = 7.5;		% pO2n		oxygen partial pressure reference
u(36) = 10;			% fvn 		vagal firing rate reference
u(37) = 5.5;		% fshn 		beta-sympth reference firing rate
u(38) = Gts;		% Gts		gain for sympathetic contribution to heart rate
u(39) = Gtv;			% Gtv		gain for vagal contribution to heart rate

% metabolic component
u(40) = c1f;		% c1f		scaling constant for fetal oxygen content
u(41) = c1f;		% c2f		scaling constant for fetal oxygen content
u(42) = alf;		% alphab 	maximum binding capacity of hemoglobin
u(43) = Hbf;		% Hbf 		fetal hemoglobin concentration
u(44) = bet;		% betad 	content of dissolved oxygen

u(45) = 9.33;		% Kf 		scaling constant for metabolic consumption - Khoo et al
u(46) = 0.068;		% cOth		threshold oxygen concentration
u(47) = 0.31;		% Omet0		initial oxygen metabolism rate

u(48) = 0.25;		% klc		kinetic constant
u(49) = 50;			% kk5		kinetic constant
u(50) = 0.04;		% kk1		kinetic constant
u(51) = 0.02;		% kk2		kinetic constant
u(52) = 0.01;		% kk4		kinetic constant
u(53) = 30;			% kk7		kinetic constant

u(54) = 0.042;		% Omet0c	initial cerebral oxygen metabolic consumption
u(55) = 0.038;		% c0thc		cerebral threshold oxygen concentration
u(56) = 1;			% KK 		KK coefficient

% additional parameters
u(57) = tauv;		% tau_v		vagal time constant in low pass filter
u(58) = taus;		% tau_s		sympathetic time constant in low pass filter
u(59) = 30;			% tauisc	
u(60) = 0.05;		% psi		constant in blood pressure normalisation integral

u(61) = 1;			% glc0		initial gluose concentration
u(62) = 0.8;		% lac0		initial lactate concentration
u(63) = 0.12;		% py0		initial pyruvate concentration
u(64) = 40;			% hc0		initial H+ ion concentration;

u(65) = 10;			% pyruvate consumption rate
u(66) = 44;			% cO2mc0	initial systemic oxygen concentration
u(67) = 0.1;		% H+ feedback to CO2
u(68) = 0.1;		% kmc00		systemic kinetic constant in metabolism

% ----------------% ----------------% ----------------% ----------------

u = u';		% transpose for matlab vector multiplicaton

% solve the system now
tic 
option = odeset('MaxStep', 1e3, 'RelTol', 1e-4);
[t, Y] = ode15s(@(t, Y)mainode1(t, Y, u), tspan, Y0, option);

toc

%% Figure plotting

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