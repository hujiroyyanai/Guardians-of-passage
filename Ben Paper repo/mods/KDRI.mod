TITLE HH k channel channel
: Hodgkin - Huxley k channel

: The model used in Melnick et al. 2003 
: implemented tau_factor as local variable in FUNCTION and used it insted of Q10 (15thMay20, Kazutaka)

NEURON {
	SUFFIX KDRI
	USEION k READ ek WRITE ik
	RANGE gkbar, ik
	GLOBAL inf
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
PARAMETER {
	v (mV)
	dt (ms)
	gkbar=0 (mho/cm2) <0,1e9>
	ek = -84 (mV)
	celsius = 37 (degC)
}
STATE {
	n h
}
ASSIGNED {
	ik (mA/cm2)
	inf[2]
	
}
LOCAL	fac[2]

INITIAL {
	rate(v*1(/mV))
	n = inf[0]
	h = inf[1]
}

BREAKPOINT {
	SOLVE states
	ik = gkbar*n*n*n*n*h*(v - ek)
}

PROCEDURE states() {	: exact when v held constant
	rate(v*1(/mV))
	n = n + fac[0]*(inf[0] - n)
	h = h + fac[1]*(inf[1] - h)
	VERBATIM
	return 0;
	ENDVERBATIM
}

UNITSOFF
FUNCTION alp(v(mV),i) { LOCAL a,b,c,tau_factor :rest=-70 order n,h
	v = v
:	q10 = 3^((celsius - 23)/10)
	tau_factor = 1
	if (i==0) {
		alp = tau_factor * 1 * .035*expM1(-v - 15, 9)
	}else if (i==1){
		alp = tau_factor * 0.0083*(1/(exp((v + 20)/10)+1) + 1)
	}
}

FUNCTION bet(v,i) { LOCAL a,b,c,tau_factor :rest=-70 order n,h
	v = v
:	q10 = 3^((celsius - 23)/10)
	tau_factor = 1
	if (i==0) {
		bet = tau_factor * 1 * .014*exp((-v + 12)/46)
	}else if (i==1){
		bet = tau_factor * 0.0083 /(exp((-v - 20)/10)+1)
	}
}

FUNCTION expM1(x,y) {
        if (fabs(x/y) < 1e-6) {
                expM1 = y*(1 - x/y/2)
        }else{
                expM1 = x/(exp(x/y) - 1)
        }
}


PROCEDURE rate(v) {LOCAL a, b, tau :rest = -70
	TABLE inf, fac DEPEND dt, celsius FROM -160 TO 100 WITH 200
	FROM i=0 TO 1 {
		a=alp(v,i)  b=bet(v,i)
		tau = 1/(a+b)
		inf[i] = a/(a+b)
		fac[i] = (1 - exp(-dt/tau))
	}
}
UNITSON
