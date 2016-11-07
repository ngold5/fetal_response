function dy = mainode1(t,y,uu)
global Hbm c1m c2m Vut Bput Rut
global T0 fv0 fvh tuco tucor ttr
%ylag1 = Z(:,1);
%ylag2 = Z(:,2);
% fetal circulation trial
% with oxygen distribution
% with heart rate included and UCO
dy = zeros(30,1);  
% model cerebral blood flow control and oxygen consumption
%{
y(1) = p_a
y(2) = p_v
y(3) = V_l

y(4) = pum
y(5) = pc

y(4+2) = pcO_a    aterial
y(5+2) = pcO_um    umbilical
y(6+2) = pcO_mc   systemic
y(7+2) = pcO_c   cerebral
 y(8+2)  regulation of Delta Tv
y(11)  = pO  mother
y(12) = regulation of Delta Ts
y(13) = fab
y(14) = fac
y(15) = fsh0
y(16) = Rcmc
y(17)  MAP
-- mother
y(18)  Pa     
y(19)  Pv
y(20)  Vl
y(21)  Pivs
--
y(22)  CO2_a
--
y(23)  CO2_um
--
y(24)  CO2_mc
--
y(25)  CO2_c
--
y(26)  fsp0
-----
y(27)  Glu
-----
y(28)  Lac
-----
y(29)  Pyr
-----
y(30)  H+
%}

% parameters
factor = 0.76/101.325;  % translate kPa ms/ml to mmHg s/ml ... 
%now pressure is in mmHg and time is in second
factor2 = 1000*factor;  % translate kPa to mmHg

factor = 1/133.322;
factor2 = 1000/133.322;

Ra = uu(1);   %6*factor;  % kPa ms/ml = 2e-5 mmHg s/cm^3
Rv = uu(2);   %30*factor*1.;
%Rum = 700*factor;
Rmc0 = uu(3);   %200*factor;
%Rc0 = 1400*factor;

Rummc = uu(4);   %698*factor;
Rumv = uu(5);   %2*factor;
Rcmc0 = uu(6);   %1398*factor;
Rcv = uu(7);   %2*factor;

Ca0 = uu(8);


Gts = uu(38);  
Gtv = uu(39);
tauv = uu(57);
taus = uu(58);

% compute pc
np = 5.9717;
mp = 4.0917;
np = 5;
mp = 4;

tmin = 0.0801;
tmax = 0.2653;
pmin = 1.0994;
pmax = 1.2028;
%pmax = 1.5;
eta = 19.0681;
phi = 1.3296;
alfa = 0.000;
theta = 1.1364;
nu = 9.1682;
%eta = 2.;
%nu = 2.;

% regulation
fevmin = uu(12);   %3.2;
%fevmin = 0.5;
fevmax = uu(13);   %6.3;
Wcv = uu(14);   %0.2;
Wbv = uu(15);    %1;
kev = uu(16);   %7.06;
fvhmin = 0;
fvhmax = 2.5;
fv01 = - 0.68;
%fv01 = -5;
kvh = uu(17);  %1.7;

fabmin = 2.52;
fabmax = 47.78;
facmin = 1.16;
facmax = 17.07;
kab = uu(18);   %4.09;
kac = uu(19);   %7.24;
%fabn = 40;
fabn = uu(20);   %25;
facn = uu(21);   %11;

%fesmax = 60;
%fesinf = 2.1;
%fes0 = 16.11;
kes = uu(22);   %0.0675;
Wby = uu(23);   %1.;
Wcy = uu(24);   %1;
%Wpy = 0;  % no pulmonary here

% add Rmc
fesmax = uu(25);    %60;
fesinf = uu(26);     %2.1;
fes0 = uu(27);   %16.11;

%fsp = fesinf + (fes0 - fesinf)*exp(kes*( -ylag2(13) + 5*ylag2(14) - ylag2(26) ));
fsp = fesinf + (fes0 - fesinf)*exp(kes*( -y(13) + 5*y(14) - y(26) ));
if(fsp>=fesmax)
    fsp = fesmax;
end

fsh = fesinf + (fes0 - fesinf)*exp(kes*( -Wby*y(13) + Wcy*y(14) - y(15) ));
%fsh = fesinf + (fes0 - fesinf)*exp(kes*( -Wby*ylag2(13) + Wcy*ylag2(14) - ylag2(15) ));

if(fsh>=fesmax)
    fsh = fesmax;
end

%do we need delay to bring back the FHR
fesmin = uu(28);   %2.66;
fs1 = 1;

if(fsh<fesmin)
    qtmp = 0;
else
    qtmp = 1;
end


a = 0.0003;
b = 6.3041;
dm = 1.1264;
cm = 1.5;
d = 5;
c0 = 0.9;  % contractility
%H = 1.3;   % heart rate

%H = 1/(T0-0.2+y(10));  % 1/s
%H = 1/(T0+y(12)+y(10));  % 1/s
%H = 1/(T0-0.177+y(12)+y(10));  % 1/s
H = 1/(T0-0.13*log(5.5-uu(28)+1)+y(12)+y(10));  % 1/s
%H = 1/(T0-0.13*log(fsh-2.66+1)+y(12)+y(10));  % 1/s
%H = 1.5;

tp = tmin + theta^nu/(H^nu+theta^nu)*(tmax-tmin);
betaH = (np+mp)/np*tp;
Pp = pmin + H^eta/(H^eta+phi^eta)*(pmax-pmin);
tt = mod(t,1/H);

f = 0; gt = 0;
if(tt<alfa)
    f = 0; gt = 0;
elseif(tt>=alfa && tt<=betaH)
    f = Pp*(tt-alfa)^np*(betaH-tt)^mp/np^np/mp^mp/((betaH-alfa)/(mp+np))^(mp+np);
elseif(tt>betaH && tt<= 1/H)
    f = 0; gt = 0;
end

%ttp = mod(tp,1/H);
ttp = tp;
fftp =0;
if(ttp<alfa)
    fftp = 0;
elseif(ttp>=alfa && ttp<=betaH)
    fftp = Pp*(ttp-alfa)^np*(betaH-ttp)^mp/np^np/mp^mp/((betaH-alfa)/(mp+np))^(mp+np);
elseif(ttp>betaH)
    fftp = 0;
end
if(abs(fftp)>0.000001)
    gt = f/fftp;
end

% add dc
dc = 0;
if(t>ttr)
dc = 2*abs(10.2*log(fsh - 2.66+1) - 8.7)/75*qtmp;
%dc = (0 + 0.5*exp((fsh-4.)/1.))/(1+exp((fsp-4)/1));
%dc = (0 + 0.5*exp(-(y(17)-pn)/1.))/(1+exp(-(y(17)-pn)/1));
end
c = c0 + dc;
%if(t>ttr+3600+1800)
%c = c0 + dc*2;
%c = 1.2;
%end

pc = a*(y(3)-b)^2 + (c*y(3)-d)*gt;

% -- mother heart
Hm = .75;
tp = tmin + theta^nu/(Hm^nu+theta^nu)*(tmax-tmin);
betaH = (np+mp)/np*tp;
Pp = pmin + Hm^eta/(Hm^eta+phi^eta)*(pmax-pmin);
tt = mod(t,1/Hm);

f = 0; gtm = 0;
if(tt<alfa)
    f = 0; gtm = 0;
elseif(tt>=alfa && tt<=betaH)
    f = Pp*(tt-alfa)^np*(betaH-tt)^mp/np^np/mp^mp/((betaH-alfa)/(mp+np))^(mp+np);
elseif(tt>betaH && tt<= 1/Hm)
    f = 0; gtm = 0;
end

%ttp = mod(tp,1/H);
ttp = tp;
fftp =0;
if(ttp<alfa)
    fftp = 0;
elseif(ttp>=alfa && ttp<=betaH)
    fftp = Pp*(ttp-alfa)^np*(betaH-ttp)^mp/np^np/mp^mp/((betaH-alfa)/(mp+np))^(mp+np);
elseif(ttp>betaH)
    fftp = 0;
end
if(abs(fftp)>0.000001)
    gtm = f/fftp;
end

pcm = (cm*y(20)-dm)*gtm;

% mother heart pressure done


% if no contraction i.e. put is not included
%Rc = Rc0;
Rmc = Rmc0;
Rcmc = Rcmc0;


quma = (y(1)-y(4))/Rummc;  % umbilical in artery
qumv = (y(4)-y(2))/Rumv;

%cOin = (cOum*qum + cOmc*qmc + cOc*qc)/qav;   % in flow


c1f = uu(40);
c2f = uu(41);
alf = uu(42);
Hbf = uu(43);
bet = uu(44);

% arterial
Sa = 1- c1f/(c1f + y(6)^3+c2f*y(6));
Sum = 1- c1f/(c1f + y(7)^3+c2f*y(7));
Smc = 1- c1f/(c1f + y(8)^3+c2f*y(8));
Sc = 1- c1f/(c1f + y(9)^3+c2f*y(9));
cOa = alf*Hbf*Sa + bet*y(6);  %  coeff 100 omitted
cOum = alf*Hbf*Sum + bet*y(7);
cOmc = alf*Hbf*Smc + bet*y(8);
cOc = alf*Hbf*Sc + bet*y(9);
%cOa = y(6);
%cOmc = y(8);


%Ca0 = 2/factor2*4;  % ml/kPa  % control amplitude of pressure

Ca = Ca0;

ptmp=1;
Ratmp = Ra;
fvh1 = fvh;
Gtv1 = Gtv;
Gts1 = Gts;
%if(ptmp==0)
%    Gts1 = Gts*(1 - fv/fv0/2);
%end
%pO20 = 14.8;   % mmHg
pO20 = uu(29);   %7.;

%if(t>ttr+3600 && t<ttr+3600*2)
%    tucor = 2*60;
%elseif(t>ttr+3600*2)
%    tucor = 60;
%end

if(t>ttr)
    if(t - floor(t/(tuco+tucor))*(tuco+tucor))<= tuco
%        if(t<ttr+(tuco+tucor)*10)
%            ptmp = 0.5;
%        elseif(t>ttr+(tuco+tucor)*10 && t<ttr+(tuco+tucor)*20)
%            ptmp = 0.25;
%        else
%            ptmp = 0;
%        end
        if(t<ttr+60*60)
            ptmp = 0.5;
        elseif(t>ttr+3600 && t<ttr+3600*2)
            ptmp = 0.25;
        elseif(t>ttr+3600*2 && t<ttr+3600*4)
           ptmp = 0;
        else
            ptmp = 1;
        end
    end
    quma = quma*ptmp;
    qumv = qumv*ptmp;
end
%    Ramin = Ra/10;
%    Ramax = Ra*10;
%   Ratmp = (Ramin + Ramax*exp(-(cOa-0.05)/0.005))/(1+exp(-(cOa-0.05)/0.005));

kcmc = uu(30);   %0.018;
cOan= uu(31);  
taucmc = 10;   
   Rcmcmax = Rcmc0;
   Rcmcmin = Rcmc0/10;
Rcmcs = (Rcmcmax + Rcmcmin*exp(-(cOa-cOan)/kcmc))/(1+exp(-(cOa-cOan)/kcmc));

dy(16) = (Rcmcs-y(16))/taucmc;

qa = (pc - y(1))/Ratmp;
qv = (y(2)-pc)/Rv;


Rmcmax = 1500*factor;
Rmcmin = 0*factor;
fspn = uu(32);  % 9
kkmc = uu(33);  % 1.44

dRmc = (Rmcmin + Rmcmax*exp((fsp-fspn)/kkmc))/(1+exp((fsp-fspn)/kkmc));
% add Rmc
Rmc = Rmc0 + dRmc;

qmc = (y(1)-y(2))/Rmc;  % systemic
%qc = (y(1)-y(2))/Rc;   % cerebral
%qca = (y(1)-y(5))/Rcmc;
qca = (y(1)-y(5))/y(16);
qcv = (y(5)-y(2))/Rcv;

qav = quma + qmc + qca; % in 3 compartments
qmv = qumv + qmc + qcv;  % out 3 compartments

% add dCv
dCv = (-10.0*log(fsp - 2.66+1)/15+0.8)*qtmp;
%Cv = 35/factor2;
%Cv = 35/factor2 + 0.5*dCv;
Cv = uu(9);   %35/factor2;
Cum = uu(10);  %11/factor2*1;
Cc = uu(11);   %0.57/factor2*1;

%qav = (pl - pav)/Rav;
%qmv = (pmv - pl)/Rmv;

% Vl
sv = 0;
sa = 0;

if(y(2)>pc)
    sv = 1;
elseif(pc>y(1))
    sa = 1;
elseif(y(2)>pc && pc>y(1))
    sa = 1; sv = 1;
end

dy(3) = sv*qv - sa*qa;  % now combined ventricle

if(sa==0)
    qa=0;
end
if(sv==0)
    qv=0;
end
    

dy(1) = (qa-qav)/Ca;
dy(2) = (qmv-qv)/Cv;

dy(4) = (quma - qumv)/Cum;
dy(5) = (qca - qcv)/Cc;

% =================== mother
Ram = 7*factor;  % kPa ms/ml = 2e-5 mmHg s/cm^3 (7,30)
Rvm = 3*factor;
Cam = 42/factor2;
Cvm = 840/factor2;
Cut = 0.8/factor2;
%Cut = 10./factor2;

Rivs = 450*factor;
Rutv = 90*factor;
Rmmc = 60*factor;

qam = (pcm - y(18))/Ram;
qvm = (y(19) - pcm)/Rvm;

qmmc = (y(18)-y(19))/Rmmc;
qutsa = (y(18)-y(21))/Rivs;
qutsv = (y(21)-y(19))/Rutv;

qavm = qmmc + qutsa;
qmvm = qmmc + qutsv;


% Vl
svm = 0;
sam = 0;

if(y(19)>pcm)
    svm = 1;
elseif(pcm>y(18))
    sam = 1;
elseif(y(19)>pcm && pcm>y(18))
    sam = 1; svm = 1;
end

if(sam==0)
    qam=0;
end
if(svm==0)
    qvm=0;
end
dy(18) = (qam - qavm)/Cam; % Pa
dy(19) = (qmvm - qvm)/Cvm;  % Pv

dy(20) = svm*qvm - sam*qam; 

dy(21) = (qutsa - qutsv)/Cut;
% ===================


Cp = [Ca;Cv;Cum;Cc];
vtmp = sum(dy(1:2).*Cp(1:2))+...
    sum(dy(4:5).*Cp(3:4)) + dy(3);
Cp2 = [Cam;Cvm;Cut];
vtmp2 = sum(dy(18:19).*Cp2(1:2))+...
    + dy(21)*Cut + dy(20);
if(abs(vtmp)>0.00001)
    disp('not conserved!-1')
elseif(abs(vtmp2)>0.00001)
    disp('not conserved!-2')
end

% add oxygen distribution, now the variables are partial pressure
%alf = 1.34e-6;
%bet = 0.031;
%Hbf = 1.7e5;
%Hbm = 1.2e5;
%c1f = 1.04e4;
%c2f = 150;
%c1m = 2.34e4;
%c2m = 150;


%cOa = y(6);
%cOmc = y(8);

Va0 = 20;    % ml
Va = Va0 + Ca*y(1);

Omet = 0;
cOin = (cOum*quma + cOmc*qmc + cOc*qca)/qav;   % in flow

dy(6) = (qa*(cOin-cOa) - Omet - cOa*Ca*dy(1))/Va...
    /(alf*Hbf*c1f*(3*y(6)^2+c2f)/(c1f+y(6)^3+c2f*y(6)).^2+bet);
%dy(6) = (qa*(cOin-cOa) - Omet - cOa*Ca*dy(1))/Va;

% ---
cCO2in = 0.56;  % mmHg  
kco2 = 0.244;
Kco2 = 0.007;
pCO2a = (y(22) - kco2)/Kco2;  % linear relation used by Khoo et al
% empirical formula
Ip = 35.5;   % mmHg
Gp = 26.5/60;  % (s*mmHg)^-1
Ve = Gp*exp(-0.05*y(6))*(pCO2a - Ip);
OMa = - Ve*y(22) + 0.1;
dy(22) = (qa*(cCO2in - y(22)) + OMa - y(22)*Ca*dy(1))/Va;  % CO2_a
%---

glc0 = uu(61);
lac0 = uu(62);
py0 = uu(63);
hc0 = uu(64);

cCO2in = (y(23)*quma + y(24)*qmc + y(25)*qca)/qav;   % in flow
GLCa = (glc0*quma + y(27)*qmc + glc0*qca)/qav;   % in flow
Laca = (lac0*quma + y(28)*qmc + lac0*qca)/qav; 
Pyra = (py0*quma + y(29)*qmc + py0*qca)/qav; 
Hc0 = (hc0*quma + y(30)*qmc + hc0*qca)/qav; 

%GLCa = glc0;   % in flow
%Laca = lac0; 
%Pyra = py0; 
%Hc0 = hc0; 


% um
Vum0 = 90;
Vum = Vum0 + Cum*y(4);
%Vum = 116;  % ml  from Sa Couto et al 2002
%Vum = Vum0;

D = 2.83e-2; % cm^3Q2/s/mmHg
Omet = 0;
Odiff = D*(y(11) - y(7));  % coupling between mother and fetus !!!
%Odiff = D*(30. - y(7));
%Odiff = 0;
%cOin_um =  0.03;  % in flow, assuming constant from mother
%cOin_um = cOum*qum/qav;   % in flow
dy(7) = (ptmp*quma*(cOin-cOum) + Odiff - Omet -...
    cOum*Cum*dy(4))/Vum...
    /(alf*Hbf*c1f*(3*y(7)^2+c2f)/(c1f+y(7)^3+c2f*y(7)).^2+bet);

%-------------------------  umbilical cord metaboic pathways
% y(23) CO2_um
%cCO2in = 0.52; %mmHg
kco2 = 0.244;
Kco2 = 0.007;
pCO2um = (y(23) - kco2)/Kco2;  % linear relation used by Khoo et al
%OMum = -D*(pCO2um - 40);
%OMum = -Odiff*0.8;
OMum = 0.1*(44-pCO2um) + 0.;
dy(23) = (ptmp*quma*(cCO2in - y(23)) + OMum - y(23)*Cum*dy(4))/Vum;  % CO2_um
%-------------------------  umbilical cord metaboic pathways

% mc
Vmc0 = 100;
Cmc = 0.85/factor2;   % ? not given in fetal paper
Cmc = 0;
%Vmc = Vmc0 + Cmc*0.5*(y(1)+y(2));
%Vmc = 220;
Vmc = Vmc0;

Omet0 = uu(47);  %3.1e-1;   % cm^3O2/s
Kf = uu(45);   %9.33;   % ml blood/s
cOth = uu(46);   %0.068;
%cOth = 0.025;
if(cOmc>= cOth)
    Omet = Omet0;
else
    Omet = Omet0 + Kf*(cOmc -cOth);
end
%cOin_mc = (cOin - cOin_um)*qmc/(qmc+qc);   % in flow, not justified
%cOin_mc = cOmc*qmc/qav;   % in flow
dy(8) = (qmc*(cOin - cOmc) - Omet - ...
    cOmc*Cmc*0.5*(dy(1)+dy(2)))/Vmc...
    /(alf*Hbf*c1f*(3*y(8)^2+c2f)/(c1f+y(8)^3+c2f*y(8)).^2+bet);
%  oxygen in unit volume
%dy(8) = (qmc*(cOin - cOmc) - Omet - ...
%    cOmc*Cmc*0.5*(dy(1)+dy(2)))/Vmc;

%-------------------------  systemic metaboic pathways
% y(24)  CO2_mc
%
kmc1 = uu(48); % klc
kmc2 = uu(49);  % kk5
kmc3 = uu(50);  % kk1
kmc4 = uu(51);  % kk2
kmc5 = uu(52);  % kk4
kmc6 = uu(53);  % kk7
kmc0 = uu(65);
kmc00 = uu(68);

khtO = uu(67);

%cCO2in = 0.56;  % mmHg  
OMmc = 0.2*Omet + khtO*(y(30)-40)^2/((y(30)-40)^2+1);  % 0.1 is due to buffer? units?
dy(24) = (qmc*(cCO2in - y(24)) + OMmc- y(24)*Cmc*0.5*(dy(1)+dy(2)))/Vmc;  % CO2_mc
pCO2mc = (y(24) - kco2)/Kco2;
%-------------------------  systemic metaboic pathways
% glucose  (mM)
dy(27) = (qmc*(GLCa - y(27)) +kmc0*y(29) - kmc00*y(27)/(y(27) + 9.) ...
    - y(27)*Cmc*0.5*(dy(1)+dy(2)))/Vmc;
% lactate  (mM)
klc = kmc1;   %0.5;
kk5 = kmc2;   %60;
%if(t>ttr+3600)
%klc = 0.6;
%kk5 = 50;
%end
dy(28) = (qmc*(Laca - y(28))+ kmc0*y(29) + klc*y(30) - y(28)*Cmc*0.5*(dy(1)+dy(2)))/Vmc;
% pyruvate  (mM)
kk1=  kmc3;  %0.04; 
kk2 = kmc4;  %0.02; 
%kk3=0.06; 
kk4 = kmc5;  %0.01;
%dy(29) = qmc/Vmc*(0.12 - y(29)) + kk1*y(27) - kk2*y(29)*y(8)/(0.12*0.1+y(29)*y(8)) ...
%    - kk3*y(29) + kk4*y(28);
dy(29) = (qmc*(Pyra - y(29)) + kk1*y(27) - kk2*cOmc - kk4*y(28) ...
    - y(29)*Cmc*0.5*(dy(1)+dy(2)))/Vmc;
% H+ (nM)
%kk5 = 35; 
%kk6 = 0.00; 
kk7 = kmc6;  %30.; 
%kk8 = 0.0005;   % notice the diff in units
dy(30) = (qmc*(Hc0 - y(30)) + kk5*y(28) + kk7*y(24) ...
    - y(30)*Cmc*0.5*(dy(1)+dy(2)))/Vmc;  % 0.1 is buffer

%-------------------------  metabolism

% c
Vc0 = 7.2;
Vc = Vc0 + Cc*y(5);
%Vc = Vc0;

Omet0c = uu(54);  %4.2e-2;
cOthc = uu(55);  %0.039;
%cOthc = 0.005;
if(cOc>= cOthc)
    Ometc = Omet0c;
else
    Ometc = Omet0c + Kf*(cOc -cOthc);
end
%cOin_c = (cOin - cOin_um)*qc/(qmc+qc);   % in flow, not justified
%cOin_c = cOc*qc/qav;   % in flow
dy(9) = (qca*(cOin - cOc) - Ometc - ...
    cOc*Cc*dy(5))/Vc...
    /(alf*Hbf*c1f*(3*y(9)^2+c2f)/(c1f+y(9)^3+c2f*y(9)).^2+bet);
%dy(9) = (qca*(cOin - cOc) - Ometc - cOc*Cc*dy(5))/Vc;

%-------------------------  cerebral compartment metaboic pathways
%cCO2in = 0.56;  % mmHg  
OMc = 0.2*Ometc;  % + .0*(y(34)-40)^2/((y(34)-40)^2+1);
dy(25) = (qca*(cCO2in - y(25)) + OMc- y(25)*Cc*dy(5))/Vc;  % CO2_c
pCO2c = (y(25) - kco2)/Kco2;  % linear relation used by Khoo et al
%---
%-------------------------  systemic metaboic pathways
% glucose  (mM)
%dy(31) = (qca*(1.2 - y(31)) - 0.1*y(31)/(y(31) + 4.) ...
%    - y(31)*Cc*dy(5))/Vc;
% lactate  (mM)
%dy(32) = (qca*(1.4 - y(32))+ 10.*y(33) + 0.1*y(34) - y(32)*Cc*dy(5))/Vc;
% pyruvate  (mM)
%kk1=0.04; kk2 = 0.02; kk3=0.06; kk4=0.01;
%dy(33) = (qca*(0.12 - y(33)) + kk1*y(31) - kk2*cOc - kk4*y(32)- y(33)*Cc*dy(5))/Vc;
% H+ (nM)
%kk5 = 40; kk6 = 0.00; kk7 = 20; kk8 = 0.0005;   % notice the diff in units
%dy(34) = (qca*(40 - y(34)) + kk5*y(32) + kk7*y(25)- y(34)*Cc*dy(5))/Vc;  %
%0.1 is buffer
%-------------------------  metabolism

%-------------------------  cerebral compartment



Sut = 1- c1m/(c1m + y(11)^3+c2m*y(11));
cOut = alf*Hbm*Sut + bet*y(11);

%qut = Bput/Rut;
qut = qutsa;
cOinut = 0.159;
cOinut = 0.2;
%Odiff = D*(y(11)-y(7));
dy(11) = (qut*(cOinut - cOut) - Odiff)/Vut...
    /(alf*Hbm*c1m*(3*y(11)^2+c2m)/(c1m+y(11)^3+c2m*y(11)).^2+bet);
%dy(11)=0;


% adding dynamic block
%tauz = 6.37;
%taup = 2.076;
%dy(17) = (y(1) + tauz*dy(1) - y(17))/taup;

% regulation on sympathetic activity
Gab = (y(13)-fabn)/(y(17)-10);
%kab = (fabmax-fabmin)/4/Gab;
Gac = (y(14)-facn)/(y(6)-pO20);
%kac = (facmax-facmin)/4/Gac;

pn = uu(34);   %40;
fabs = (fabmin + fabmax*exp((y(17)-pn)/kab))/(1+exp((y(17)-pn)/kab));
%fabs = (fabmax + fabmin*exp((y(6)-pO20)/kab))/(1+exp((y(6)-pO20)/kab));
%y(13) = fabs;
%--CO2 modification
%KK = uu(56)*log(pCO2mc/44.0)+1;
KK = uu(56)*log(pCO2mc/uu(66))+1;
%--
facs = KK*(facmax + facmin*exp((y(6)-pO20)/kac))/(1+exp((y(6)-pO20)/kac));

%facs = (facmax + facmin*exp((y(6)-pO20)/kac))/(1+exp((y(6)-pO20)/kac));

%fabs = (fabmin + fabmax*exp((y(13)-fabn)/kab))/(1+exp((y(13)-fabn)/kab));
%facs = (facmin + facmax*exp((y(14)-facn)/kac))/(1+exp((y(14)-facn)/kac));

%fv = fv0 + fvh*(1-y(9)/pO20); % no delay 
% modified model
fvh1 = (fvhmin + fvhmax*exp(-(y(9)-pO20)/kvh))/(1+exp(-(y(9)-pO20)/kvh)) - fv01;

fv = Wbv*(fevmin + fevmax*exp((y(13)-fabn)/kev))/(1+exp((y(13)-fabn)/kev)) ...
    + Wcv*y(14) + fvh1;
%fvh1 = (fvhmin + fvhmax*exp(-(ylag1(9)-pO20)/kvh))/(1+exp(-(ylag1(9)-pO20)/kvh)) - fv01;
%fv = Wbv*(fevmin + fevmax*exp((ylag1(13)-fabn)/kev))/(1+exp((ylag1(13)-fabn)/kev)) ...
%    + Wcv*ylag1(14) + fvh1;
%fv = fv0 + fvh*(1-Z(9,1)/pO20);  


fvn = uu(36);
fshn = uu(37);

Gtv1 = Gtv + (Gtv/10*exp((fsh-fshn)/.5))./(1+exp((fsh-fshn)/.5));
dy(10) = (Gtv1*fv - y(10))/tauv;  % zero delay... 


%if(ptmp==0)
%  Gts1 = Gts/2.d0
%end

%Gts1 = (Gts/10+Gts*exp(-(fv-10)/0.05))./(1+exp(-(fv-10)/0.05));
%Gts1 = (Gts/20+0.16 + (Gts+0.16)*exp(-(fv-10)/0.01))./(1+exp(-(fv-10)/0.01));
Gts1 = (Gts+0.13 + (Gts/20+0.13)*exp((fv-fvn)/0.05))./(1+exp((fv-fvn)/0.05));
if(Gts1<0.01)
Gts1=0.01;
end


DtTs = (Gts1*log((fsh - fesmin + fs1)/fs1)+...
    0.5*Gts*log((fsh - fesmin + fs1)/(5.571-fesmin+fs1)))*qtmp;       
%DtTs = y(16)*log((fsh - fesmin + fs1)/fs1)*qtmp;     
dy(12) = (DtTs - y(12))/taus;


tauab = 2;
tauac = 1.5;
% fab  y(13)
dy(13) = (fabs - y(13))/tauab;
% fac  y(14)
dy(14) = (facs - y(14))/tauac;
% fsh,0  y(15)  % threshold

tauisc = uu(59);  %30;

fsh0min = -49.38;
fsh0max = 3.59;
fsh0n = 11.36;
ksh0 = 6;
pO2n = uu(35);   %7.5;

fshs = (fsh0min + fsh0max*exp((y(6)-pO2n)/ksh0))/(1+exp((y(6)-pO2n)/ksh0));
fsp0min = 7.33;
fsp0max = 13.32;
fsp0n = 7.52;
ksp0 = 2;
fsps = (fsp0min + fsp0max*exp((y(6)-pO2n)/ksp0))/(1+exp((y(6)-pO2n)/ksp0));

dy(15) = (fshs - y(15))/tauisc;
dy(26) = (fsps - y(26))/tauisc;

% logistic equation for Gts
%taur = 1.;
%dy(16) = (y(16)-Gts1)*(Gts*(1-fv/fv0/2)-y(16))*taur;
%dy(16) = 0;

% equation for average firing rate
psi = uu(60);  %0.05;
N = (1-exp(-psi*t))/psi;
dy(17) = (y(1) - y(17))/N;

%-------------------


return