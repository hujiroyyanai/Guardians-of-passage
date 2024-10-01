/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
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
 
#define nrn_init _nrn_init__B_Na
#define _nrn_initial _nrn_initial__B_Na
#define nrn_cur _nrn_cur__B_Na
#define _nrn_current _nrn_current__B_Na
#define nrn_jacob _nrn_jacob__B_Na
#define nrn_state _nrn_state__B_Na
#define _net_receive _net_receive__B_Na 
#define _f_rates _f_rates__B_Na 
#define rates rates__B_Na 
#define states states__B_Na 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define alpha_shift _p[0]
#define alpha_shift_columnindex 0
#define beta_shift _p[1]
#define beta_shift_columnindex 1
#define tau_factor _p[2]
#define tau_factor_columnindex 2
#define ina _p[3]
#define ina_columnindex 3
#define gnabar _p[4]
#define gnabar_columnindex 4
#define tadj _p[5]
#define tadj_columnindex 5
#define inf (_p + 6)
#define inf_columnindex 6
#define tau (_p + 8)
#define tau_columnindex 8
#define m _p[10]
#define m_columnindex 10
#define h _p[11]
#define h_columnindex 11
#define ena _p[12]
#define ena_columnindex 12
#define Dm _p[13]
#define Dm_columnindex 13
#define Dh _p[14]
#define Dh_columnindex 14
#define a (_p + 15)
#define a_columnindex 15
#define b (_p + 17)
#define b_columnindex 17
#define v _p[19]
#define v_columnindex 19
#define _g _p[20]
#define _g_columnindex 20
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
 
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
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_alpha(void);
 static void _hoc_beta(void);
 static void _hoc_rates(void);
 static void _hoc_trap(void);
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
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_B_Na", _hoc_setdata,
 "alpha_B_Na", _hoc_alpha,
 "beta_B_Na", _hoc_beta,
 "rates_B_Na", _hoc_rates,
 "trap_B_Na", _hoc_trap,
 0, 0
};
#define alpha alpha_B_Na
#define beta beta_B_Na
#define trap trap_B_Na
 extern double alpha( _threadargsprotocomma_ double , double );
 extern double beta( _threadargsprotocomma_ double , double );
 extern double trap( _threadargsprotocomma_ double , double );
 
static void _check_rates(double*, Datum*, Datum*, NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, int _type) {
   _check_rates(_p, _ppvar, _thread, _nt);
 }
 /* declare global and static user variables */
#define usetable usetable_B_Na
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_B_Na", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "alpha_shift_B_Na", "mV",
 "beta_shift_B_Na", "mV",
 "ina_B_Na", "mA/cm2",
 "gnabar_B_Na", "mho/cm2",
 "tadj_B_Na", "1",
 "inf_B_Na", "1",
 "tau_B_Na", "ms",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "usetable_B_Na", &usetable_B_Na,
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
"B_Na",
 "alpha_shift_B_Na",
 "beta_shift_B_Na",
 "tau_factor_B_Na",
 0,
 "ina_B_Na",
 "gnabar_B_Na",
 "tadj_B_Na",
 "inf_B_Na[2]",
 "tau_B_Na[2]",
 0,
 "m_B_Na",
 "h_B_Na",
 0,
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 21, _prop);
 	/*initialize range parameters*/
 	alpha_shift = 0;
 	beta_shift = 0;
 	tau_factor = 1;
 	_prop->param = _p;
 	_prop->param_size = 21;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
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

 void _B_NA_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 21, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 B_Na /Users/royyanai/Documents/repos/repo_with_michael/mods/B_NA.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_inf[2];
 static double *_t_tau[2];
static int _reset;
static char *modelname = "HH sodium channel";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_rates(_threadargsprotocomma_ double);
static int rates(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_rates(_threadargsprotocomma_ double _lv);
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargscomma_ v ) ;
   Dm = ( inf [ 0 ] - m ) / tau [ 0 ] ;
   Dh = ( inf [ 1 ] - h ) / tau [ 1 ] ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau[0] )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau[1] )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
   rates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau[0])))*(- ( ( ( inf[0] ) ) / tau[0] ) / ( ( ( ( - 1.0 ) ) ) / tau[0] ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau[1])))*(- ( ( ( inf[1] ) ) / tau[1] ) / ( ( ( ( - 1.0 ) ) ) / tau[1] ) - h) ;
   }
  return 0;
}
 
double alpha ( _threadargsprotocomma_ double _lv , double _li ) {
   double _lalpha;
 if ( _li  == 0.0 ) {
     _lalpha = 0.182 * trap ( _threadargscomma_ - _lv + 7.0 - 35.0 + alpha_shift , 9.0 ) ;
     }
   else if ( _li  == 1.0 ) {
     _lalpha = 0.061 * trap ( _threadargscomma_ - _lv + 13.0 - 48.0 + alpha_shift , 3.0 ) + 0.0166 ;
     }
   
return _lalpha;
 }
 
static void _hoc_alpha(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  alpha ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double beta ( _threadargsprotocomma_ double _lv , double _li ) {
   double _lbeta;
 if ( _li  == 0.0 ) {
     _lbeta = 0.124 * trap ( _threadargscomma_ _lv - 7.0 + 35.0 + beta_shift , 9.0 ) ;
     }
   else if ( _li  == 1.0 ) {
     _lbeta = 0.0018 * trap ( _threadargscomma_ _lv - 13.0 + 84.0 + beta_shift , 18.0 ) ;
     }
   
return _lbeta;
 }
 
static void _hoc_beta(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  beta ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double trap ( _threadargsprotocomma_ double _lx , double _ly ) {
   double _ltrap;
 if ( fabs ( _lx / _ly ) < 1e-6 ) {
     _ltrap = _ly * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _ltrap = _lx / ( exp ( _lx / _ly ) - 1.0 ) ;
     }
   
return _ltrap;
 }
 
static void _hoc_trap(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  trap ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 static double _mfac_rates, _tmin_rates;
  static void _check_rates(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_tadj;
  if (!usetable) {return;}
  if (_sav_tadj != tadj) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rates =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_rates)/200.; _mfac_rates = 1./_dx;
   for (_i=0, _x=_tmin_rates; _i < 201; _x += _dx, _i++) {
    _f_rates(_p, _ppvar, _thread, _nt, _x);
    for (_j = 0; _j < 2; _j++) { _t_inf[_j][_i] = inf[_j];
}    for (_j = 0; _j < 2; _j++) { _t_tau[_j][_i] = tau[_j];
}   }
   _sav_tadj = tadj;
  }
 }

 static int rates(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _lv) { 
#if 0
_check_rates(_p, _ppvar, _thread, _nt);
#endif
 _n_rates(_p, _ppvar, _thread, _nt, _lv);
 return 0;
 }

 static void _n_rates(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rates(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_rates * (_lv - _tmin_rates);
 if (isnan(_xi)) {
  for (_j = 0; _j < 2; _j++) { inf[_j] = _xi;
}  for (_j = 0; _j < 2; _j++) { tau[_j] = _xi;
}  return;
 }
 if (_xi <= 0.) {
 for (_j = 0; _j < 2; _j++) { inf[_j] = _t_inf[_j][0];
} for (_j = 0; _j < 2; _j++) { tau[_j] = _t_tau[_j][0];
} return; }
 if (_xi >= 200.) {
 for (_j = 0; _j < 2; _j++) { inf[_j] = _t_inf[_j][200];
} for (_j = 0; _j < 2; _j++) { tau[_j] = _t_tau[_j][200];
} return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 for (_j = 0; _j < 2; _j++) {double *_t = _t_inf[_j]; inf[_j] = _t[_i] + _theta*(_t[_i+1] - _t[_i]);}
 for (_j = 0; _j < 2; _j++) {double *_t = _t_tau[_j]; tau[_j] = _t[_i] + _theta*(_t[_i+1] - _t[_i]);}
 }

 
static int  _f_rates ( _threadargsprotocomma_ double _lv ) {
   {int  _li ;for ( _li = 0 ; _li <= 1 ; _li ++ ) {
     a [ _li ] = alpha ( _threadargscomma_ _lv , ((double) _li ) ) ;
     b [ _li ] = beta ( _threadargscomma_ _lv , ((double) _li ) ) ;
     tau [ _li ] = 1.0 / ( a [ _li ] + b [ _li ] ) / tau_factor ;
     } }
   inf [ 0 ] = a [ 0 ] / ( a [ 0 ] + b [ 0 ] ) ;
   inf [ 1 ] = 1.0 / ( 1.0 + exp ( ( _lv + 75.0 - 11.0 ) / 9.0 ) ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_rates(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  m = m0;
 {
   tadj = pow( 3.0 , ( ( celsius - 23.0 ) / 10.0 ) ) ;
   rates ( _threadargscomma_ v ) ;
   m = inf [ 0 ] ;
   h = inf [ 1 ] ;
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];

#if 0
 _check_rates(_p, _ppvar, _thread, _nt);
#endif
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
  ena = _ion_ena;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ina = gnabar * m * m * m * h * ( v - ena ) ;
   }
 _current += ina;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  ena = _ion_ena;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  ena = _ion_ena;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = m_columnindex;  _dlist1[0] = Dm_columnindex;
 _slist1[1] = h_columnindex;  _dlist1[1] = Dh_columnindex;
  for (_i=0; _i < 2; _i++) {  _t_inf[_i] = makevector(201*sizeof(double)); }
  for (_i=0; _i < 2; _i++) {  _t_tau[_i] = makevector(201*sizeof(double)); }
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/royyanai/Documents/repos/repo_with_michael/mods/B_NA.mod";
static const char* nmodl_file_text = 
  "TITLE HH sodium channel\n"
  ": Hodgkin - Huxley squid sodium channel\n"
  "\n"
  ": The model used in Melnick et al. 2004 Adapt 5 and 11 mV\n"
  ": moved alpha/beta_shift from RANGE to GLOBAL (15thMay20, Kazutaka)\n"
  ": implemented tau_factor to RANGE and used it insted of tadj (15thMay20, Kazutaka)\n"
  ": (should set alpha/beta_shift as GLOBAL variable when use test_EXinitial)\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX B_Na\n"
  "	USEION na READ ena WRITE ina\n"
  "	RANGE gnabar, ina\n"
  "	RANGE inf, tau\n"
  "	RANGE tadj, tau_factor, alpha_shift, beta_shift\n"
  "	: GLOBAL alpha_shift, beta_shift\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "}\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "PARAMETER {\n"
  "	ena = 53 (mV)\n"
  "	alpha_shift = 0 (mV)\n"
  "	beta_shift = 0 (mV)\n"
  "	tau_factor = 1\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	m h\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	celsius (degC)\n"
  "	v (mV)\n"
  "	ina (mA/cm2)\n"
  "\n"
  "	gnabar (mho/cm2)\n"
  "	tadj (1)\n"
  "\n"
  "	inf[2] (1)\n"
  "	tau[2] (ms)\n"
  "	a[2] (1/ms)\n"
  "	b[2] (1/ms)	\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	tadj = 3^((celsius - 23) / 10)\n"
  "	rates(v)\n"
  "	m = inf[0]\n"
  "	h = inf[1]\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD cnexp\n"
  "	ina = gnabar*m*m*m*h*(v - ena)\n"
  "}\n"
  "\n"
  "DERIVATIVE states{\n"
  "	rates(v)\n"
  "	m' = (inf[0] - m) / tau[0]\n"
  "	h' = (inf[1] - h) / tau[1]\n"
  "}\n"
  "\n"
  "FUNCTION alpha(v(mV),i) {\n"
  "	if       (i==0){\n"
  "		alpha = 0.182*trap(-v + 7 - 35 + alpha_shift, 9)\n"
  "	}else if (i==1){\n"
  "		alpha = 0.061*trap(-v + 13 - 48 + alpha_shift, 3) + 0.0166\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION beta(v,i) {\n"
  "	if       (i==0){\n"
  "		beta = 0.124 * trap(v - 7 + 35 + beta_shift, 9)\n"
  "	}else if (i==1){\n"
  "		beta = 0.0018 * trap(v - 13 + 84 + beta_shift, 18)\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION trap(x,y) {\n"
  "	if (fabs(x/y) < 1e-6) {\n"
  "		trap = y*(1 - x/y/2)\n"
  "	}else{\n"
  "		trap = x/(exp(x/y) - 1)\n"
  "	}\n"
  "}\n"
  "\n"
  "PROCEDURE rates(v) {\n"
  "	TABLE inf, tau DEPEND tadj FROM -100 TO 100 WITH 200\n"
  "	FROM i=0 TO 1 {\n"
  "		a[i] = alpha(v,i) b[i]=beta(v,i)\n"
  "		tau[i] = 1 / (a[i] + b[i]) / tau_factor\n"
  "	}\n"
  "	inf[0] = a[0] / (a[0] + b[0])\n"
  "	inf[1] = 1 / (1 + exp((v + 75 - 11) / 9))\n"
  "}\n"
  ;
#endif
