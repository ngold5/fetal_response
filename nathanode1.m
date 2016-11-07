function dy = nathanode1(t, y, uu)
global Vut Bput Rut
global T0 fv0 fvh tuco tuco ttr
% Nathan's reduced version of Qiming's model
% fetal circulation - hemodynamics
% oxygen distribution and CO2 distribution
% heart rate included and UCO
% regulation model
% original metabolic model included
dy = zeros(16,1)
%{
	y(1) = umbilical pressure p_U
	y(2) = fetal arterial pressure p_A
	y(3) = fetal venous pressure p_V
	y(4) = heart volume change V_H

	Metabolics
	y(5) = pcO_A	arterial partial pressure of O2
	y(6) = pcO_U	umbilical partial pressure of O2
	y(7) = pcO_F	fetal partial pressure of O2	
	y(8) = CO2_U	umbilical CO2 concentration
	y(9) = CO2_F	fetal carbon dioxide concentration

%}

% parameters
factor = 1 / 133.322;
factor2 = 1000 * factor;

Ra = uu(1);		% 6*factor, arterial resistnace
Rv = uu(2);		% 30*factor, venous resistance
Rmc0 = uu(3);	% 200*factor, systemic fetal resitance
Rummc = uu(4);	% 698*factor
Rumv = u(5);	% 2*factor, umbilical venous resitance
Rcmc0 = uu(6);	% 1398*factor, cerebral resistance ??
Rcv = uu(7);	% 2*factor, cerebebral venous resitance ??

Ca0 = uu(8);	% arterial capacitance

% regulation parameters
Gts = uu(38);	% sympathetic gain
Gtv = uu(39);	% vagal gain
tauv = uu(57);	% vagal time constant
taus = uu(58);	% sympathetic time constant

% heart volume change model
np = 5;
mp = 4;

tmin = 0.0801;
tmax = 0.2653;
pmin = 1.0994;
pmax = 1.2028;

eta = 19.0681;
phi = 1.3296;
alfa = 0;
theta = 1.1364;
chi = 9.1682;

% regulation
% to be properly filled in


% heart rate properties
a = 0.0003;
b = 6.3041;
dm = 1.1264;
cm = 1.5;		% maternal cardiac contractility
d = 5;
c0 = 0.9;		% baseline fetal heart contractility

% calculate heart rate H
H = 1 / (T0 - 0.13*log(5.5 - uu(28) + 1) + y(12) + y(10)); % 1/s
% y(10) accounts for the vagal contribution to heart rate
% y(12) accounts for the sympathetic contribution to heart rate
% these need to be changed in our model and simplified

% calculate constituent functions for heart as a pump
tp = tmin + theta^chi/(H^chi + theta^chi) * (tmax - tmin);
betaH = (np + mp)/np * tp;
Pp = pmin + H^eta/(H^eta + phi^eta) * (pmax - pmin);
tt = mod(t,1/H);

% activation function
f = 0; gt = 0;
if (tt < alfa)
	f = 0; gt = 0;
elseif (tt >= alfa && tt <= betaH)
	f = Pp*(tt - alfa)^np * (betaH - tt)^mp/np^np/mp^mp/((betaH - alfa)/(mp+np))^(mp+np);
elseif (tt > betaH && tt <= 1/H)
	f = 0; gt = 0;
end

% activation function, fH(tp) where ttp as input - denominator
ttp = tp;
fttp = 0;
if (ttp < alfa)
	fttp = 0;
elseif (ttp >= alfa && ttp <= betaH)
	fttp = Pp*(ttp - alfa) * (betaH - ttp)^mp/np^np/mp^mp/((betaH-alfa)/(mp+np))^(mp+np);
elseif (ttp > betaH)
	fttp = 0;
end

if(abs(fttp) > 0.0000001)	% prevent division by zero singularity
	gt = f / fttp;	% scaled activation function for input to fetal heart
end

% --------------- CHECK THIS -----------------
% variation in fetal heart contractility - still to be determined properly
% currently placehold using Qiming's input values - need modify from
% beta sympathetic activity to just sympathetic, should be ok to do
dc = 0; % initial offset parameter
%if (t > ttr) % we are in the contractility transition state
%	dc = 2 * abs(10.2 * log(fsh - 2.66 + 1) - 8.7) / 75 * qtmp;
%end

% Change dynamics for contractility to piecewise functions dependent
% on UCO - this can be modified to incorporate sympathetic firing 
% information, alternatively

if (t > ttr)
	if ((t - floor(t/(tuco + tucor)) * (tuco + tucor)) <= tuco)
		if (t < ttr + 3600)
			dc = 0.0256;
		elseif (t > ttr + 3600 && t < ttr + 3600*2)
			dc = 0.0638;
		elseif (t > ttr + 3600*2 && t < ttr + 3600*4)
			dc = 0.652;
		end
	end

c = c0 + dc;
% --------------- CHECK THIS -----------------

%heart pressure-volume relationship
pc = a*(y(4) - b)^2 + (c*y(4) - d)*gt;

% --------------- CHECK THIS -----------------

% --------------------------------------------------------------
% hemodynamic flow equations
Rf = Rmc0; % set initial systemic resistance
Rcmc = Rcmc0; % set initial cerebral systemic resitance

% umbilical flows:
qua = (y(2) - y(1)) / Rummc; % umbilical arterial flow qUA
quv = (y(4) - y(3)) / Rumv;	% umbilical venous flow qUV

% oxygen concentration parameters
c1f = uu(40);
c2f = uu(41);
alf = uu(42);
Hbf = uu(43);
bet = uu(44);

% arterial conversions to pressures and concentrations
% define denominators for relation between oxygen concentration and
% oxygen partial pressure
Sa = 1 - c1f/(c1f + y(5)^3 + c2f*y(5));	% arterial denominator
Su = 1 - c1f/(c1f + y(6)^3 + c2f*y(6));	% umbilical denominator
Sf = 1 - c1f/(c1f + y(7)^3 + c2f*y(7)); % fetal systemic denominator

% oxygen concentration concentration conversion from partial pressure
cOa = alf*Hbf*Sa + bet*y(5); % arterial O2 concentration
cOu = alf*Hbf*Su + bet*y(6); % umbilical O2 concentration
cOf = alf*Hbf*Sf + bet*y(7); % fetal systemic O2 concentration

Ca = Ca0;	% initial value for arterial capicitance

% UCO parameters and set up now
ptmp = 1;	% initial flow multiplier for umbilical pressure
Ratmp = Ra;	% arterial resistance
fvh1 = fvh;	% initial vagal firing rate
Gtv1 = Gtv;	% vagal gain parameter for FHR
Gts1 = Gts;	% sympathetic gain parameter for FHR

% initalise pressure
pO20 = uu(29);	% mmHg pressure variable

% ----------------------------------------------------------------
% UCO

% ttr = time to recovery
% tuco = time of UCO - duration of occlusion
% tucor = time of UCO recovery - duration of recovery between occlusion

if (t > ttr)
	if ((t - floor(t/(tuco + tucor)) * (tuco + tucor) <= tuco) % in UCO now
		% set up UCO lengths and modify flow allowance from UC
		if (t < ttr + 3600) % first hour of UCO mild
			ptmp = 0.5; % mild occlusion - only reduce flow to half
		elseif (t > ttr + 3600 && t < ttr + 3600*2) % second UCO - moderate
			ptmp = 0.25;	% 1/4 of allowed flow - moderate UCO
		elseif (t > ttr + 3600*2 && t < ttr + 3600*4) % third UCO - severe
		  	ptmp = 0; % no allowed flow during this time - severe UCO series
		else % no UCO
			ptmp = 1; % flow is normal
		end
	end
	qua = qua * ptmp; % umbilical artery flow
	quv = quv * ptmp; % umbilical venous flow
end
			
qf = (y(2) - y(3)) / Rf;	% fetal systemic flow
qa = qua + qf; % arterial flow - inflow of umbilical and fetal compartment
qv = quv + qf; % venous flow - outflow of umbilical and fetal compartment

% add compliances - analog of capicitance
Cv = uu(9);	% venous compliance
Cum = uu(10); % umbilical compliance - should this vary???

% heart pumping action for flow - indicator function
sv = 0;	% venous flow - inflow to heart
sa = 0; % arterial flow - outflow of heart

if (y(3) > pc) % pv higher than heart pressure
	sv = 1; % inflow to heart
elseif (pc > y(2)) % pa higher than heart pressure
	sa = 1; % outflow from heart
elseif (y(3) > pc && pc > y(1)) % combined inflow outflow
	sa = 1; sv = 1;
end

dy(4) = sv*qv - sa*qa;	% combined ventricle volume change

if (sa == 0)
	qa = 0;	% no outflow
end
if (sv ==0)
	qv = 0; %  no inflow
end

% hemodynamic equations
dy(2) = (qa - )


			







