% Qiming's ODE model
% Qiming's model in his code uses the same variable names as Van der Hout-
% van der Jagt
% This will translate them into their variable names as they are actually
% used in our paper
% ------------------------------------------------------------------------
% Fetal circulation
% with oxygen distribution, heart rate included, and UCO
% ------------------------------------------------------------------------

function dy = qimingode(t, y, uu)

% declare global variables
global Hbm c1m c2m Vut Bput Rut
global T0 fv0 fvh tuco tucor ttr

% ylag1 = Z(:, 1); 		Qiming's model omits the time delay lag terms
% ylag2 = Z(:, 2);		Qiming's model omits the time delay lag terms

% initialise memory for ODE system
dy = zeros(30, 1);

% variables
%{
	y(1) = p_a		arterial blood pressure
	y(2) = p_v		venous blood pressure
	y(3) = V_l		compartmental volume
	y(4) = pum		umbilical pressure
	y(5) = pc		cerebral blood pressure
	y(6) = pcO_a	arterial oxygen pressure
	y(7) = pcO_um	umbilical oxygen pressure
	y(8) = pcO_mc	systemic oxygen pressure
	y(9) = pcO_O 	cerebral oxygen pressure
	y(10) = regulation of Delta Tv
	y(11) = pO mother
	y(12) = regulation of Delta Ts
	y(13) = fab		afferent firing rate baroreceptor
	y(14) = fac		afferent firing rate chemoreceptor
	y(15) = fsh0	threshold firing rate for sympathetic activity to heart
	y(16) = Rcmc	cerebral flow resistance
	y(17) = MAP		mean absolute blood pressure
	--------------
	Mother
	--------------
	y(18) = Pa 		maternal arterial blood pressure
	y(19) = Pv 		maternal venous blood pressure
	y(20) = Vl 		maternal compartmental volume
	y(21) = Pivs	intervillous space pressure
	--------------
	Metabolics
	--------------
	y(22) = CO2_a	arterial oxygen concentration
	y(23) = CO2_um	umbilical oxygen concentration
	y(24) = CO2_mc	system oxygen concentration
	y(25) = CO2_c 	cerebral oxygen concentration
	-------------
	y(26) = fsp0	threshold firing rate for sympathetic activity to vessels
	--------------
	Metabolics
	--------------
	y(27) = Glu		glucose concentration
	y(28) = Lac		lactate concentration
	y(29) = Pyr		pyruvate concentration
	y(30) = H+		hydrogen ion concentration
%}
% ----------------------------------------------------------------
% parameters
% ----------------------------------------------------------------
factor = 0.76 / 101.325; % translate kPa ms/ml to mmHg s/ml
% this converts pressure to mmHg and time into seconds
factor2 = 1000 * factor; % translate kPa to mmHg

factor = 1/133.322;
factor2 = 1000 / 133.322;

Ra = uu(1);		% 6 * factor; % kPa ms/ml = 2e-5 mmHg s/cm^3; arterial resistance
Rv = uu(2); 	% 30 * factor*1.; venous resistance
Rmc0 = uu(3); 	% 200 * factor; systemic resistance

Rummc = uu(4);	% 698 * factor; umbilical cord resistance
Rumv = uu(5);	% 2 * factor; umbilical venous resistance
Rcmc0 = uu(6);	% 1398 * factor; cerebral threshold resistance
Rcv = uu(7);	% cerebral venous resistance

Ca0 = uu(8);

Gts = uu(38);	% Ts gain - sympathetic contribution to heart rate gain
Gtv = uu(39);	% Tv gain - vagal (parasympth) contribution to heart rate gain
tauv = uu(57);	% time constant for vagal activity - parasympathetic
taus = uu(58);	% time constant for sympathetic activity

% ----------------
% compute pc
% constants in cardiovascular model
% ----------------
np = 5.9717; 
mp = 4.0917;
np = 5;
mp = 4;

tmin = 0.0801; 
tmax = 0.2653;
pmin = 1.0994;
pmax = 1.2028;

eta = 19.0681;
phi = 1.3296;
alfa = 0.0003;
theta = 1.1364;
nu = 9.1682;

% ----------------
% regulation
% ----------------
% parameters for vagal nerve firing
fevmin = uu(12);	% 3.2; vagal minimum firing rate 
fevmax = uu(13);	% 6.3; vagal maximum firing rate
Wcv = uu(14);		% 0.2; chemoreceptor weighting parameter
Wbv = uu(15);		% 1; baroreceptor weighting parameter
kev = uu(16);		% 7.06; parameter to determine slope of response in sigmoid

% fvh
fvhmin = 0;
fvhmax = 2.5;
fv01 = -0.68;
kvh = uu(17);	% 1.7; parameter to determine slope of response in sigmoid

% baro and chemoreceptor parameters
fabmin = 2.52;	% afferent firing rate - baroreceptor - minimum
fabmax = 47.78; % afferent firing rate - baroreceptor - maximum
facmin = 1.16; 	% afferent firing rate - chemoreceptor - minimum
facmax = 17.07;	% afferent firing rate - chemoreceptor = minimum
kab = uu(18);	% 4.09; control parameter in baroreceptor sigmoid
kac = uu(19);	% 7.24; control parameter in chemoreceptor sigmoid

% sympathetic nervous parameters
kes = uu(22);		% 0.0675;
Wby = uu(23);		% 1; baroreceptor weighting factor
Wcy = uu(24);		% 1; chemoreceptor weighting factor
fesmax = uu(25);	% 60; 
fesinf = uu(26);	% 2.1;
fes0 = uu(27);		% 16.11;

% Calculate firing rate for sympathetic innervation now
% fsp = alpha sympathetic peripheral firing rate
fsp = fesinf + (fes0 - fesinf) * exp(kes * (-y(13) + 5*y(14) - y(26)));
if (fsp >= fspmax)
	fsp = fspmax; % limit the maximum sypathetic firing rate
end

% fsh = beta sympathetic heart firing rate
fsh = fesinf + (fes0 - fesinf) * exp(kes * (-Wby*y(13) + Wcy*y(14) - y(15)));
if (fsh >= fshmax)
	fsh = fshmax;
end

fesmin = uu(28);	% 2.66
fs1 = 1;

if (fsh < fesmin)
	qtmp = 0;
else
	qtmp = 1;
end

% Heart rate H
% Heart rate is given as a frequency per minute
% So we must convert the heart period T into the frequency
% we keep the heart rate in 1/s for computational convenience - it is the same
H = (T0 - 0.13*log(5.5 - uu(28) + 1) + y(12) + y(10))^(-1); % 1/s
% Equation 18 in van der Hout-van der Jagt


% Now consider the pressure-volume relationship in the fetal heart
% This is taken from Olufsen
% Regulation equations that feed into cardiovascular model - from Olufsen
tp = tmin + theta^nu / (H^nu + theta^nu) * (tmax - tmin);
betaH = (np + mp)/np * tp;
Pp = pmin + H^eta /(H^eta + phi^eta) * (pmax - pmin);
tt = mod(t, 1/H); % make sure time converts properly
% gh(t) = fh(t) / fh(tp)
f = 0; gt = 0;
if (tt < alfa) % indicator
	f = 0; gt = 0;
elseif (tt > alfa && tt <= betaH)
	f = Pp * (tt - alfa)^np * (betaH - ttp)^mp/np^np/mp^mp/((betaH - alfa)/(mp+np))^(mp+np);
elseif (ttp > betaH)
	fftp = 0;
end
if (abs(fftp) > 0.000001)
	gt = f / fftp;
end


% ----------------
% add dc
% ----------------
dc = 0;
if (t > ttr)
	dc = 2 * abs(10.2 * log(fsh - 2.66 + 1) - 8.7) / 75*qtmp;
end
c = c0 + dc;
% What exactly is dc doing??? This is a very strange part of the code
% ----------------

% equation (5) from Qiming's model
pc = a*(y(3) - b)^2 + (c*y(3) - d)*gt;

% ----------------
% mother's heart now
% ----------------
Hm = 0.75;
tp = tmin + theta^nu / (Hm^nu + theta^nu) * (tmax - tmin);
betaH = (np + mp)/np * tp;
Pp = pmin + Hm^eta/(Hm^eta + phi^eta)*(pmax - pmin);
tt = mod(t, 1/Hm);

f = 0; gtm = 0;
if (tt < alfa)
	f = 0; gtm = 0;
elseif (tt >= alfa && tt <= betaH)
	f = Pp*(tt-alfa)^np*(betaH-tt)^mp/np^np/mp^mp/((betaH-alfa)/(mp+np))^(mp+np);
elseif (tt>betaH && tt <= 1/Hm)
	f = 0; gtm = 0;
end

ttp = tp;
fftp = 0;
if (ttp < alfa)
	fftp = 0;
elseif (ttp >= alfa && ttp <= betaH)
	fftp = Pp*(ttp - alfa)^np*(betaH - ttp)^mp/np^np/mp^np/((betaH - alfa)/(mp + np))^(mp + np);
elseif (ttp> betaH)
	fftp = 0;
end

if (abs(fttp) > 0.000001)
	gtm = f / fttp;
end

pcm = (cm * y(20) - dm) * gtm;

% conclusion of mother heart pressure

% ----------------
% if no contraction, put - uterine pressure - is not included

Rmc = Rmc0;		% systemic resistance
Rcmc = Rcmc0;	% cerebral resistance

quma = (y(1) - y(4)) / Rummc;	% artery blood flow in umbilical cord
qumv = (y(4) - y(2)) / Rumv;	% venous blood flow in umbilical cord

% ----------------
% oxygen concentration parameters
% ----------------
c1f = uu(40);	% fetal scaling constant
c2f = uu(41);	% fetal scaling constant
a1f = uu(42);	% maximum binding capacity of hemoglobin
Hbf = uu(43);	% fetal hemoglobin concetration
bet = uu(44); 	% beta constant content of disolved oxygen

% ----------------
% arterial 
% these equations below are for oxygen diffusion
% they are modified equations from Couto et al 2002
% ----------------
Sa = 1 - c1f / (c1f + y(6)^3 + c2f*y(6));		% arterial
Sum = 1 - c1f / (c1f + y(7)^3 + c2f*y(7));		% umbilical
Smc = 1 - c1f / (c1f + y(8)^3 + c2f*y(8));		% systemic
Sc = 1 - c1f / (c1f + y(9)^3 + c2f*y(9));		% cerebral
% ----------------
% equations below are from Couto for total oxygen content in blood
% ----------------
cOa = alf*Hbf*Sa + bet*y(6); % equation 15 from Couto - arterial
cOum = alf*Hbf*Sum + bet*y(7); % equation 15 Couto - umbilical
cOmc = alf*Hbf*Smc + bet*y(8); % equation 15 Couto - systemic
cOc = alf*Hbf*Sc + bet*y(9); % equation 15 Couto - cerebral


Ca = Ca0;

ptmp = 1;
Ratmp = Ra;	
fvh1 = fvh;		% vagal firing rate
Gtsv1 = Gtv;	% vagal gain
Gts1 = Gts;		% sympathetic gain

% -----------------------------------------------
% Occlusions UCO
% -----------------------------------------------
if (t > ttr)
	if (t - floor(t/(tuco + tucor)) * (tuco + tucor)) <= tuco
		if (t < ttr + 60*60) % 1 min UCO
			ptmp = 0.5;
		elseif (t>ttr + 3600 && t < ttr + 3600*2) % 2 min UCO
			ptmp = 0.25;
		elseif (t > ttr + 3600*2 && t < ttr + 3600 * 4) % open
			ptmp = 1;
		end
	end
	quma = quma*ptmp;
	qumv = qumv*ptmp;
end
	
% autoregulation component
kcmc = uu(30);	% 0.018; steepness of cerebral vascular resistance sigmoidal curve
cOan = uu(31);
taucmc = 10; 	% time constant for autoregulation
Rcmcmax = Rcmc0;
Rcmcmin = Rcmc0 / 10;
% define the reference intermediate sigmoidal curve for the low pass filter ODE
Rcmcs = (Rcmcmax + Rcmcmin*exp(-(cOa - cOan)/kcmc)) / (1 + exp(-(cOa - cOan)/kcmc));

% define low pass filter - first order ODE
% cerebral resistance low pass filter
dy(16) = (Rcmcs - y(16)) / taucmc;

qa = (pc - y(1)) / Ratmp;
qv = (y(2) - pc) / Rv;

% system autoregulation component
Rmcmax = 1500*factor;
Rmcmin = 0*factor;
fspn = uu(32);	% 9
kkmc = uu(33);	% 1.44

dRmc = (Rmcmin + Rmcmax * exp((fsp - fspn)/kkmc)) / (1 + exp((fsp - fspn) / kkmc));

Rmc = Rmc0 + dRmc;

qmc = (y(1) - y(2)) / Rmc;
qca = (y(1) - y(5)) / y(16);
qcv = (y(5) - y(2)) / Rcv;

% Blood flow branching - mass conservations - equations (6) (7) Wang
qav = quma + qmc + qca; % in 3 compartments
qmv = qumv + qmc + qcv;	% out 3 compartments

% ----------------
% venous compliance
% ----------------
dCv = (-10.0 * log(fsp - 1.66)/15.8) * qtmp;

Cv = uu(9);
Cum = uu(10);
Cc = uu(11);

% Vl
sv = 0; sa = 0;
% indicator functions for heart pumping
if (y(2) > pc)
	sv = 1;
elseif(pc > y(1))
	sa = 1;
elseif(y(2) > pc && pc > y(1))
	sa = 1; sv = 1;
end

dy(3) = sv*qv - sa*qa; % combined ventricle - equation (4) in Wang
% the above parts count as indicator functions that are scanned through
% the indicators decide whether or not the heart will pump blood
% or of the ventricle is closed as to receive blood

% indicator flow functions
if (sa == 0)
	qa = 0; % no arterial flow
end
if (sv ==0)
	qv = 0; % no venous flow
end

dy(1) = (qa - qav) / Ca;		% arterial blood pressure ODE
dy(2) = (qmv - qv) / Cv; 		% venous blood pressure ODE
dy(4) = (quma - qumv) / Cum;	% umbilical blood pressure ODE
dy(5) = (qca - qcv) / Cc;		% cerebral blood pressure ODE

% ----------------% ----------------% ----------------
% ----------------% ----------------% ----------------
% MOTHER
% ----------------% ----------------% ----------------
% ----------------% ----------------% ----------------
Ram = 7 * factor;
Rvm = 3 * factor;
Cam = 42 / factor2;
Cvm = 840 / factor2;
Cut = 0.8 / factor2;

Rivs = 450 * factor;
Rutv = 90 * factor;
Rmmc = 60 * factor;

qam = (pcm - y(18)) / Ram;			% maternal arterial flow
qvm = (y(19) - pcm) / Rvm;			% maternal venous flow

qmmc = (y(18) - y(19)) / Rmmc;		% maternal systemic flow
qutsa = (y(18) - y(21)) / Rivs;		% maternal uterine arterial flow
qutsv = (y(21) - y(19)) / Rutsv;	% maternal uterine venous flow

% impose mass conservation on flow into uterine compartment
qavm = qmmc + qutsa;
qmvm = qmmc + qutsv;


% Vl - volume
svm = 0; sam = 0;

% flow indicator functions for uterine and placental fow
if (y(19) > pcm)
	svm = 1;
elseif(pcm > y(18))
	sam = 1;
elseif(y(19) > pcm && pcm > y(18))
	svm = 1; sam = 1;
end

% logic control flow for maternal flows
if (sam == 0)
	qam = 0;
end
if (svm == 0)
	qvm = 0;
end

% maternal pressure ODEs
dy(18) = (qam - qavm) / Cam;	% Pa - maternal arterial pressure
dy(19) = (qmvm - qvm) / Cvm;	% Pv - maternal venous pressure

% heart blood volume ODe
dy(20) = svm * qvm - sam * qam;

% uterine pressure ODE
dy(21) = (qutsa - qutsv) / Cut;

% ----------------% ----------------% ----------------% ----------------
% Assert conservation of total volume
% Gather compliances - fetal
Cp = [Ca; Cv; Cum; Cc];
% volume conservation
vtmp = sum(dy(1:2).*Cp(1:2)) + sum(dy(4:5).*Cp(3:4)) + dy(3);

% gather compliances - maternal
Cp2 = [Cam; Cvm; Cut];
% assert total volume conseration
vtmp2 = sum(dy(18:19).*Cp2(1:2)) + dy(21)*Cut + dy(20);

if (abs(vtmp) > 0.00001)
	disp('not conserved!-1')
elseif (abs(vtmp2)>0.00001)
	disp('not conserved!-2')
end


% ----------------% ----------------% ----------------% ----------------
% ----------------% ----------------% ----------------% ----------------
% Metabolic component
% ----------------% ----------------% ----------------% ----------------
% ----------------% ----------------% ----------------% ----------------

% O2
Va0 = 20;
Va = Va0 + y(1);

Omet = 0;		% O2 metabolic rates
cOin = (cOum*quma + cOmc*qmc + cOc*qca) / qav; % feeding oxygen concentration into the fetus

% oxygen concentration in compartments
dy(6) = (qa * (cOin - cOa) - Omet - cOa*Ca*dy(1)) / Va ...
	/ (alf*Hbf*c1f*(3*y(6)^2 + c2f)/(c1f+y(6)^3+c2f*y(6)).^2 + bet);

% ----------------% ----------------% ----------------% ----------------	
% CO2
cCO2in = 0.56;	% mmHg
kco2 = 0.244;
Kco2 = 0.007;
% implement Khoo et al linear relationship between O2 and CO2 eq (23)
pCO2a = (y(22) - kco2) / Kco2;	% emperical relationship

Ip = 35.5;		% mmHg
Gp = 26.5 / 60; % (s*mmHg)^(-1)
Ve = Gp*exp(-0.05*y(6)) * (pCO2a - Ip);
OMa = -Ve*y(22) + 0.1;	% metabolic production rate

% arterial CO2_a concentration
dy(22) = (qa * (cCO2in - y(22)) + OMa - y(22)*Ca*dy(1)) / Va;

% ----------------% ----------------% ----------------% ----------------
% Metabolites
% Glucose (glc); Pyruvate (py); Lactate (lac); H+ ions (hc)

% intial values of metabolites
glc0 = uu(61);		% glucose 
lac0 = uu(62);		% lactate
py0 = uu(63);		% pyruvate
hc0 = uu(64);		% H+ ions

% metabolite concentrations
cCO2in = (y(23) * quma + y(24) * qmc + y(25) * qca) / qav; % inflow of CO2
GlCa = (glc0*quma + y(27)*qumc + glc0*qca) / qav; % inflow of glucose
Laca = (lac0*quma + y(28)*qumc + lac0*qca) / qav; % inflow of lactate
Pyra = (py0*quma + y(29)*qumc + py0*qca) / qav; % inflow of pyruvate
Hc0 = (hc0*quma + y(30)*qumc + hc0*qca) / qav; % inflow of H+ ions

Vum0 = 90;	% initial umbilical volume
Vum =Vum0 + Cum*y(4); 

D = 2.83e-2;	% diffusion coefficient 
Omet = 0;		% oxygen metabolic consumption in umbilical cord
Odiff = D * (y(11) - y(7)); % maternal and fetal coupling!!! 
% Odiff describes the oxygen diffusion in the placenta as determined
% by the oxygen partial pressure difference between intervillous space
% and the umbilical cord - equation (16)

dy(7) = (ptmp*quma*(cOin - cOum) + Odiff - Omet - ...
	cOum*Cum*dy(4)) / Vum ...
	/ (alf*Hbf*c1f*(3*y(7)^2+c2f) / (c1f+y(7)^3+c2f*y(7)).^2 + bet);

% ----------------% ----------------% ----------------% ----------------
% Umbilical cord metabolic pathways
% ----------------% ----------------% ----------------% ----------------
pC02um = (y(23) - kco2) / Kco2; 	% Khoo et al linear relationship

% OMum = -D*(pCO2um - 40);
% OMum = -Odiff*0.8;
OMum = 0.1*(44 - pCO2um) + 0.;
% CO2_um - umbilical cord CO2 concentration
dy(23) = (ptmp*quma*(cCO2in - y(23)) + OMum - y(23)*Cum*dy(4)) / Vum;
% ----------------% ----------------% ----------------% ----------------

% ----------------% ----------------% ----------------% ----------------
% mc - systemic components
% ----------------% ----------------% ----------------% ----------------
Vmc0 = 100;		% initial systemic volume
Cmc = 0.85 / factor2;	% not given in fetal paper?
Cmc = 0;
Vmc = Vmc0;		% set to initial value

Omet0 = uu(47);	% 3.1e-1; cm^302/s
Kf = uu(45);	% 9.33; ml blood/s
cOth = uu(46);	% 0.068	O2 threshold to switch metabolic consumption to
% linear function of oxygen concentration
% piecewise function for metabolic oxygen consumption
if (cOmc > cOth)
	Omet = Omet0;
else 
	Omet = Omet0 + Kf*(cOmc - cOth);
end

% systemic oxygen concentration
% oxygen in unit volume
dy(8) = (qmc*(cOin - cOmc) - Omet - ...
	cOmc*Cmc*0.5*(dy(1)+dy(2))) / Vmc ...
	/ (alf*Hbf*c1f*(3*y(8)^2*c2f)/(c1f+y(8)^3+c2f*y(8)).^2 + bet);
% ----------------% ----------------% ----------------% ----------------

% Metabolite kinetics for mass conservation ODEs

kmc1 = uu(48);	% klc
kmc2 = uu(49);	% kk5
kmc3 = uu(50);	% kk1
kmc4 = uu(51);	% kk2
kmc5 = uu(52);	% kk4
kmc6 = uu(53);	% kk7
kmc0 = uu(65);
kmc00 = uu(68);
% ----------------% ----------------% ----------------% ----------------
khtO = uu(67);

OMmc = 0.2*met + khtO*(y(30) - 40)^2 / ((y(30) - 40)^2 + 1);
% systemic metabolic pathway
dy(24) = (qmc*(cCO2in - y (24)) + OMmc - y(24)*Cmc*0.5*(dy(1)+dy(2))) / Vmc;
% CO2 systemic - pC02mc - pressure equation
pC02mc = (y(24) - kc02) / Kc02;

% ----------------% ----------------% ----------------% ----------------
% glucose
% ----------------% ----------------% ----------------% ----------------
dy(27) = (qmc*(GLCa - y(27)) + kmc0*y(29) - kmc00*y(27)/(y(27) + 9.) ...
	- y(27)*Cmc*0.5*(dy(1)+dy(2))) / Vmc;	% glucose
% ----------------% ----------------% ----------------% ----------------
% lactate
% ----------------% ----------------% ----------------% ----------------
klc = kmc1;
kk5 = kmc2;
dy(28) = (qmc*(Laca - y(28)) + kmc0*y(29) + klc*y(30) - y(28)*Cmc*0.5*(dy(1)+dy(2))) / Vmc; % lactate

% ----------------% ----------------% ----------------% ----------------
% Pyruvate
% ----------------% ----------------% ----------------% ----------------
kk1 = kmc3;
kk2 = kmc4; 
kk4 = kmc5;

dy(29) = (qmc*(Pyra - y(29)) + kk1*y(27) - kk2*cOmc - kk4*y(28) ...
	- y(29)*Cmc*0.5*(dy(1) + dy(2))) / Vmc;		% pyruvate

% ----------------% ----------------% ----------------% ----------------
% H+ ions
% ----------------% ----------------% ----------------% ----------------
kk7 = kmc6;

dy(30) = (qmc*(Hc0 - y(30)) + kk5*y(28) + kk7*y(24) ...
	- y(30)*Cmc*0.5*(dy(1) + dy(2))) / Vmc;

% ----------------% ---------------- 

Vc0 = 7.2;
Vc = Vc0 + Cc*y(5);
if (cOc >= cOthc)
	Ometc = Ometc0;
else
	Ometc = Ometc + Kf*(cOc - c0thc);
end

% cerebral oxygen concentration
dy(9) = (qca*(cOin - cOc) - Ometc - ...
	cOc*Cc*dy(5)) / Vc ...
	/ (alf*hbf*c1f*(3*y(9)^2 + c2f) / (c1f+y(9)^3+c2f*y(9)).^2 + bet);

% ----------------% ----------------% ----------------% ----------------
% ---------------- cerebral metabolic pathways

OMc = 0.2*Ometc;
% cerebral CO2 concentration - cCO2_c
dy(25) = (qca*(cCO2in - y(25)) + OMc - y(25)*Cc*dy(5)) / Vc; % CO2_c
% pCO2c - cerebral CO2 concentration
pCO2c = (y(25) - kc02) / Kc02;

% ----------------% ----------------% ----------------% ----------------
% ----------------% ----------------% ----------------% ----------------
% Cerebral compartment
% ----------------% ----------------% ----------------% ----------------
% ----------------% ----------------% ----------------% ----------------
Sut = 1 - c1m / (c1m + y(11)^3 + c2m*y(11));
cOut = alf*Hbm*Sut + bet*y(11);

qut = qutsa;
cOinut = 0.159;
cOinut = 0.2;

% oxygen concentration diffusion from mother to fetus
dy(11) = (qut*(cOinut - cOut) - Odiff) / Vut...
	/ (alf*Hbm*c1m*(3*y(11)^2 + c2m) / (c1m + y(11)^2 + c2m*y(11)).^2 + bet);
% ----------------% ----------------% ----------------% ----------------
% ----------------% ----------------% ----------------% ----------------
% ----------------% ----------------% ----------------% ----------------

% regulation on sympathetic activity
Gab = (y(13) - fabn) / (y(17) - 10);	% baroreceptor
Gac = (y(14) - facn) / (y(6) - pO20);	% chemoreceptor

pn = uu(34);
% baroreceptor sympathetic firing star - sigmoid transfer function
fabs = (fabmin + fabmax*exp((y(17) - pn)/kb)) / (1+exp((y(17)-pn)/kab));

% ----------------% ----------------
% chemoreceptor firing rates now
KK = uu(56) * log(pCO2mc/uu(66)) + 1;
% chemoreceptor sympathetic firing rate star - sigmoid transfer function
fcrs = KK * (fcrmax + fcrmin*exp((y(6) - pO20)/kac)) / (1+exp((y(6) - pO20)/kac));

% ----------------% ----------------
% vagal firing rate vh1 component
fvh1 = (fvhmin + fvhmax*exp(-(y(9)-pO20)/kvh)) / (1+exp(-(y(9) - pO20)/kvh)) - fv01;
% vagal firing rate
fv = Wbv * (fevmin + fevmax*exp((y(13) - fabn)/kev)) / (1+exp((y(13) - fabn)/kev)) ...
	+ Wcv*y(14) + fvh1;

fvn = uu(36);
fshb = uu(37);

% vagal gain
Gtv1 = Gtv + (Gtv / 10*exp((fsh - fshn)/.5))./(1+exp((fsh - fshn)./5));
% vagal contribution to heart rate
dy(10) = (Gtv1*fv - y(10)) / tauv;

% sympathetic gain
Gts1 = (Gts + 0.13 + (Gts/20 + 0.13)*exp((fv - fvn)/0.05))./(1+exp((fv-fvn)/0.05));
if (Gts1 < 0.01)
	Gts1 = 0.01;
end

% delay term - this is entirely from van der Jagt-van der Hout et al
% delay term for sympathetic contribution to heart rate period
DtTs = (Gts1*log((fsh - fesmin + fs1)/fs1) + ...
	0.5*Gts*log((fsh - fesmin + fs1) / (5.571-fesmin+fs1)))*qtmp;

% sympathetic activity in heart rate
dy(12) = (DtTs - y(12)) / taus;

tauab = 2;		% baroreceptor time constant in low pass filter - ODE
tauac = 1.5;	% chemoreceptor time constant in low pass filter - ODE

% low pass filters for baro and chemoreceptors
dy(13) = (fabs - y(13)) / tauab; 	% baroreceptor filter
dy(14) = (facs - y(14)) / tauac;	% chemoreceptor filter


% ----------------% ----------------
% beta sympathetic activity
tauisc = uu(59);		% control parameter for alpha and beta sympathetic response

fsh0min = -49.38;
fsh0max = 3.59;
fsh0n = 11.36;
ksh0 = 6;
pO2n = uu(35);

% sigmoidal transfer function for beta sympathetic firing rate star
fshs = (fsh0min + fsh0max*exp((y(6) - pO2n)/ksh0)) / (1 + exp((y(6) - pO2n)/ksh0));

% alpha sympathetic activity
fsp0min = 7.33;
fsp0max = 13.32;
fsp0n = 7.52;
ksp0 = 2;

% sigmoid transfer function for alpha sympathetic firing rate star
fsps = (fsp0min + fsp0max*exp((y(6) - pO2n)/ksp0)) / (1+exp((y(6) - pO2n)/ksp0));

% low pass filters for beta and alpha sympathetic firing rates - ODE
dy(15) = (fshs - y(15)) / tauisc;		% beta sympathetic firing rate
dy(26) = (fsps - y(26)) / tauisc;		% alpha sympathetic firing rate


% equation for average firing rate
% compute the integral analytically and just evaluate numerically
psi = uu(60);
N = (1 - exp(-psi*t))/psi;
dy(17) = (y(1) - y(17)) / N;
