/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__borgka
#define _nrn_initial _nrn_initial__borgka
#define nrn_cur _nrn_cur__borgka
#define _nrn_current _nrn_current__borgka
#define nrn_jacob _nrn_jacob__borgka
#define nrn_state _nrn_state__borgka
#define _net_receive _net_receive__borgka 
#define rates rates__borgka 
#define states states__borgka 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gkabar _p[0]
#define gkabar_columnindex 0
#define vhalfn _p[1]
#define vhalfn_columnindex 1
#define vhalfl _p[2]
#define vhalfl_columnindex 2
#define vhalfm _p[3]
#define vhalfm_columnindex 3
#define vhalfk _p[4]
#define vhalfk_columnindex 4
#define a0l_fast _p[5]
#define a0l_fast_columnindex 5
#define a0l_slow _p[6]
#define a0l_slow_columnindex 6
#define a0n _p[7]
#define a0n_columnindex 7
#define fsr _p[8]
#define fsr_columnindex 8
#define zetan _p[9]
#define zetan_columnindex 9
#define zetal _p[10]
#define zetal_columnindex 10
#define gmn _p[11]
#define gmn_columnindex 11
#define gml _p[12]
#define gml_columnindex 12
#define zetam _p[13]
#define zetam_columnindex 13
#define zetak _p[14]
#define zetak_columnindex 14
#define gka _p[15]
#define gka_columnindex 15
#define n _p[16]
#define n_columnindex 16
#define l_fast _p[17]
#define l_fast_columnindex 17
#define l_slow _p[18]
#define l_slow_columnindex 18
#define ek _p[19]
#define ek_columnindex 19
#define Dn _p[20]
#define Dn_columnindex 20
#define Dl_fast _p[21]
#define Dl_fast_columnindex 21
#define Dl_slow _p[22]
#define Dl_slow_columnindex 22
#define ik _p[23]
#define ik_columnindex 23
#define l _p[24]
#define l_columnindex 24
#define _g _p[25]
#define _g_columnindex 25
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_alpm(void);
 static void _hoc_alpl(void);
 static void _hoc_alpk(void);
 static void _hoc_alpn(void);
 static void _hoc_betl(void);
 static void _hoc_betn(void);
 static void _hoc_rates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_borgka", _hoc_setdata,
 "alpm_borgka", _hoc_alpm,
 "alpl_borgka", _hoc_alpl,
 "alpk_borgka", _hoc_alpk,
 "alpn_borgka", _hoc_alpn,
 "betl_borgka", _hoc_betl,
 "betn_borgka", _hoc_betn,
 "rates_borgka", _hoc_rates,
 0, 0
};
#define alpm alpm_borgka
#define alpl alpl_borgka
#define alpk alpk_borgka
#define alpn alpn_borgka
#define betl betl_borgka
#define betn betn_borgka
 extern double alpm( double );
 extern double alpl( double );
 extern double alpk( double );
 extern double alpn( double );
 extern double betl( double );
 extern double betn( double );
 /* declare global and static user variables */
#define linf linf_borgka
 double linf = 0;
#define ninf ninf_borgka
 double ninf = 0;
#define taun taun_borgka
 double taun = 0;
#define taul_slow taul_slow_borgka
 double taul_slow = 0;
#define taul_fast taul_fast_borgka
 double taul_fast = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gkabar_borgka", "mho/cm2",
 "vhalfn_borgka", "mV",
 "vhalfl_borgka", "mV",
 "vhalfm_borgka", "mV",
 "vhalfk_borgka", "mV",
 "a0l_fast_borgka", "/ms",
 "a0l_slow_borgka", "/ms",
 "a0n_borgka", "/ms",
 "fsr_borgka", "1",
 "zetan_borgka", "1",
 "zetal_borgka", "1",
 "gmn_borgka", "1",
 "gml_borgka", "1",
 "zetam_borgka", "1",
 "zetak_borgka", "1",
 0,0
};
 static double delta_t = 0.01;
 static double l_slow0 = 0;
 static double l_fast0 = 0;
 static double n0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "ninf_borgka", &ninf_borgka,
 "linf_borgka", &linf_borgka,
 "taul_fast_borgka", &taul_fast_borgka,
 "taul_slow_borgka", &taul_slow_borgka,
 "taun_borgka", &taun_borgka,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"borgka",
 "gkabar_borgka",
 "vhalfn_borgka",
 "vhalfl_borgka",
 "vhalfm_borgka",
 "vhalfk_borgka",
 "a0l_fast_borgka",
 "a0l_slow_borgka",
 "a0n_borgka",
 "fsr_borgka",
 "zetan_borgka",
 "zetal_borgka",
 "gmn_borgka",
 "gml_borgka",
 "zetam_borgka",
 "zetak_borgka",
 0,
 "gka_borgka",
 0,
 "n_borgka",
 "l_fast_borgka",
 "l_slow_borgka",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 26, _prop);
 	/*initialize range parameters*/
 	gkabar = 0.012;
 	vhalfn = -45;
 	vhalfl = -67;
 	vhalfm = -67;
 	vhalfk = -45;
 	a0l_fast = 0.023;
 	a0l_slow = 0.023;
 	a0n = 0.04;
 	fsr = 0.6;
 	zetan = -4;
 	zetal = 2;
 	gmn = 0.45;
 	gml = 1;
 	zetam = 4;
 	zetak = -5;
 	_prop->param = _p;
 	_prop->param_size = 26;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _borgka_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 26, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 borgka /Users/royyanai/Documents/repos/repo_with_michael/mods/borgka.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Borg-Graham type generic K-A channel for a Sympathetic Preganglionic Neuron";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[3], _dlist1[3];
 static int states(_threadargsproto_);
 
double alpn (  double _lv ) {
   double _lalpn;
 _lalpn = exp ( 1.e-3 * zetan * ( _lv - vhalfn ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lalpn;
 }
 
static void _hoc_alpn(void) {
  double _r;
   _r =  alpn (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double alpk (  double _lv ) {
   double _lalpk;
 _lalpk = exp ( 1.e-3 * zetak * ( _lv - vhalfk ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lalpk;
 }
 
static void _hoc_alpk(void) {
  double _r;
   _r =  alpk (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double betn (  double _lv ) {
   double _lbetn;
 _lbetn = exp ( 1.e-3 * zetan * gmn * ( _lv - vhalfn ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lbetn;
 }
 
static void _hoc_betn(void) {
  double _r;
   _r =  betn (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double alpl (  double _lv ) {
   double _lalpl;
 _lalpl = exp ( 1.e-3 * zetal * ( _lv - vhalfl ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lalpl;
 }
 
static void _hoc_alpl(void) {
  double _r;
   _r =  alpl (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double alpm (  double _lv ) {
   double _lalpm;
 _lalpm = exp ( 1.e-3 * zetam * ( _lv - vhalfm ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lalpm;
 }
 
static void _hoc_alpm(void) {
  double _r;
   _r =  alpm (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double betl (  double _lv ) {
   double _lbetl;
 _lbetl = exp ( 1.e-3 * zetal * gml * ( _lv - vhalfl ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lbetl;
 }
 
static void _hoc_betl(void) {
  double _r;
   _r =  betl (  *getarg(1) );
 hoc_retpushx(_r);
}
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   Dn = ( ninf - n ) / taun ;
   Dl_fast = ( linf - l ) / taul_fast ;
   Dl_slow = ( linf - l ) / taul_slow ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v ) ;
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / taun )) ;
 Dl_fast = Dl_fast  / (1. - dt*( 0.0 )) ;
 Dl_slow = Dl_slow  / (1. - dt*( 0.0 )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
    n = n + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / taun)))*(- ( ( ( ninf ) ) / taun ) / ( ( ( ( - 1.0 ) ) ) / taun ) - n) ;
    l_fast = l_fast - dt*(- ( ( ( linf - l ) ) / taul_fast ) ) ;
    l_slow = l_slow - dt*(- ( ( ( linf - l ) ) / taul_slow ) ) ;
   }
  return 0;
}
 
static int  rates (  double _lv ) {
   double _la , _lq10 , _lb ;
 _lq10 = pow( 3.0 , ( ( celsius - 30.0 ) / 10.0 ) ) ;
   _la = alpn ( _threadargscomma_ _lv ) ;
   _lb = alpk ( _threadargscomma_ _lv ) ;
   ninf = 1.0 / ( 1.0 + _lb ) ;
   taun = betn ( _threadargscomma_ _lv ) / ( _lq10 * a0n * ( 1.0 + _la ) ) ;
   _la = alpl ( _threadargscomma_ _lv ) ;
   _lb = alpm ( _threadargscomma_ _lv ) ;
   linf = 1.0 / ( 1.0 + _lb ) ;
   taul_fast = betl ( _threadargscomma_ _lv ) / ( _lq10 * a0l_fast * ( 1.0 + _la ) ) ;
   taul_slow = betl ( _threadargscomma_ _lv ) / ( _lq10 * a0l_slow * ( 1.0 + _la ) ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  l_slow = l_slow0;
  l_fast = l_fast0;
  n = n0;
 {
   rates ( _threadargscomma_ v ) ;
   n = ninf ;
   l_fast = linf ;
   l_slow = linf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ek = _ion_ek;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   l = fsr * l_fast + ( 1.0 - fsr ) * l_slow ;
   gka = gkabar * n * l ;
   ik = gka * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ek = _ion_ek;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ek = _ion_ek;
 { error =  states();
 if(error){fprintf(stderr,"at line 80 in file borgka.mod:\n	SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = n_columnindex;  _dlist1[0] = Dn_columnindex;
 _slist1[1] = l_fast_columnindex;  _dlist1[1] = Dl_fast_columnindex;
 _slist1[2] = l_slow_columnindex;  _dlist1[2] = Dl_slow_columnindex;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/royyanai/Documents/repos/repo_with_michael/mods/borgka.mod";
static const char* nmodl_file_text = 
  "TITLE Borg-Graham type generic K-A channel for a Sympathetic Preganglionic Neuron\n"
  "\n"
  "COMMENT\n"
  "	Description: A-type transient K current for a Sympathetic Preganglionic Neuron.	\n"
  "	Author: Linford Briant\n"
  "	\n"
  "	A-type transient K current = \"IA\"\n"
  "	\n"
  "	Sympathetic Preganglionic Neurones = \"SPNs\"\n"
  "	\n"
  "	Note that this is a modified version of IA found widely on SenseLab. This version has had steady-state kinetics \n"
  "	that have been fit to data for the IA in SPN according to Whyment et al. (2011).\n"
  "	\n"
  "	Whyment et al. (2011), PMID: 21211550\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	v 		(mV)\n"
  "	ek 		(mV)\n"
  "	celsius 	(degC)\n"
  "	gkabar=0.012 	(mho/cm2)\n"
  "	vhalfn=-45	(mV)\n"
  "	vhalfl=-67	(mV)\n"
  "	vhalfm=-67	(mV)\n"
  "	vhalfk=-45	(mV)\n"
  "	a0l_fast=0.023	(/ms)\n"
  "	a0l_slow=0.023	(/ms)\n"
  "	a0n=0.04	(/ms)\n"
  "	fsr=0.6	(1)\n"
  "	zetan=-4	(1)\n"
  "	zetal=2    	(1)\n"
  "	gmn=0.45   	(1)\n"
  "	gml=1 	  	(1)\n"
  "	zetam=4		(1)\n"
  "	zetak=-5	(1)\n"
  "}\n"
  "\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX borgka\n"
  "	USEION k READ ek WRITE ik\n"
  "        RANGE gkabar,gka,vhalfn,vhalfl,a0l_fast,a0l_slow,a0n,zetan,zetal,gmn,gml,zetam,zetak,vhalfm,vhalfk,fsr\n"
  "        GLOBAL ninf,linf,taul_fast,taul_slow,taun\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	n\n"
  "	l_fast\n"
  "	l_slow\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "        rates(v)\n"
  "        n=ninf\n"
  "        l_fast=linf\n"
  "		l_slow=linf\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	ik (mA/cm2)\n"
  "        ninf\n"
  "        linf      \n"
  "        taul_fast\n"
  "		taul_slow\n"
  "        taun\n"
  "        gka\n"
  "		l\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD cnexp\n"
  "	l=fsr*l_fast+(1-fsr)*l_slow\n"
  "	gka = gkabar*n*l\n"
  "	ik = gka*(v-ek)\n"
  "}\n"
  "\n"
  "\n"
  "FUNCTION alpn(v(mV)) {\n"
  "  alpn = exp(1.e-3*zetan*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) \n"
  "}\n"
  "\n"
  "FUNCTION alpk(v(mV)) {\n"
  "  alpk = exp(1.e-3*zetak*(v-vhalfk)*9.648e4/(8.315*(273.16+celsius))) \n"
  "}\n"
  "\n"
  "FUNCTION betn(v(mV)) {\n"
  "  betn = exp(1.e-3*zetan*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) \n"
  "}\n"
  "\n"
  "FUNCTION alpl(v(mV)) {\n"
  "  alpl = exp(1.e-3*zetal*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))) \n"
  "}\n"
  "\n"
  "FUNCTION alpm(v(mV)) {\n"
  "  alpm = exp(1.e-3*zetam*(v-vhalfm)*9.648e4/(8.315*(273.16+celsius))) \n"
  "}\n"
  "\n"
  "FUNCTION betl(v(mV)) {\n"
  "  betl = exp(1.e-3*zetal*gml*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))) \n"
  "}\n"
  "\n"
  "DERIVATIVE states { \n"
  "        rates(v)\n"
  "        n' = (ninf - n)/taun\n"
  "        l_fast' = (linf - l)/taul_fast\n"
  "		l_slow' = (linf - l)/taul_slow\n"
  "}\n"
  "\n"
  "PROCEDURE rates(v (mV)) { :callable from hoc\n"
  "        LOCAL a,q10,b\n"
  "        q10=3^((celsius-30)/10)\n"
  "        a = alpn(v)\n"
  "	b = alpk(v)\n"
  "        ninf = 1/(1 + b)\n"
  "        taun = betn(v)/(q10*a0n*(1 + a))\n"
  "        a = alpl(v)\n"
  "	b = alpm(v)\n"
  "        linf = 1/(1 + b)\n"
  "        taul_fast= betl(v)/(q10*a0l_fast*(1 + a))\n"
  "		taul_slow= betl(v)/(q10*a0l_slow*(1 + a))\n"
  "}\n"
  ;
#endif
