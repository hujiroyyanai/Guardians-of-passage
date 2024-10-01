#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _AMPA_DynSyn_reg(void);
extern void _B_A_reg(void);
extern void _B_DR_reg(void);
extern void _B_NA_reg(void);
extern void _CaIntraCellDyn_reg(void);
extern void _GABAa_DynSyn_reg(void);
extern void _GABAb_DynSyn_reg(void);
extern void _Glycine_DynSyn_reg(void);
extern void _HH2_reg(void);
extern void _HH2new_reg(void);
extern void _KDR_reg(void);
extern void _KDRI_reg(void);
extern void _NK1_DynSyn_reg(void);
extern void _NMDA_DynSyn_reg(void);
extern void _SS_reg(void);
extern void _borgka_reg(void);
extern void _iCaAN_reg(void);
extern void _iCaL_reg(void);
extern void _iKCa_reg(void);
extern void _iNaP_reg(void);
extern void _vecevent_reg(void);
extern void _vsource_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"mods/AMPA_DynSyn.mod\"");
    fprintf(stderr, " \"mods/B_A.mod\"");
    fprintf(stderr, " \"mods/B_DR.mod\"");
    fprintf(stderr, " \"mods/B_NA.mod\"");
    fprintf(stderr, " \"mods/CaIntraCellDyn.mod\"");
    fprintf(stderr, " \"mods/GABAa_DynSyn.mod\"");
    fprintf(stderr, " \"mods/GABAb_DynSyn.mod\"");
    fprintf(stderr, " \"mods/Glycine_DynSyn.mod\"");
    fprintf(stderr, " \"mods/HH2.mod\"");
    fprintf(stderr, " \"mods/HH2new.mod\"");
    fprintf(stderr, " \"mods/KDR.mod\"");
    fprintf(stderr, " \"mods/KDRI.mod\"");
    fprintf(stderr, " \"mods/NK1_DynSyn.mod\"");
    fprintf(stderr, " \"mods/NMDA_DynSyn.mod\"");
    fprintf(stderr, " \"mods/SS.mod\"");
    fprintf(stderr, " \"mods/borgka.mod\"");
    fprintf(stderr, " \"mods/iCaAN.mod\"");
    fprintf(stderr, " \"mods/iCaL.mod\"");
    fprintf(stderr, " \"mods/iKCa.mod\"");
    fprintf(stderr, " \"mods/iNaP.mod\"");
    fprintf(stderr, " \"mods/vecevent.mod\"");
    fprintf(stderr, " \"mods/vsource.mod\"");
    fprintf(stderr, "\n");
  }
  _AMPA_DynSyn_reg();
  _B_A_reg();
  _B_DR_reg();
  _B_NA_reg();
  _CaIntraCellDyn_reg();
  _GABAa_DynSyn_reg();
  _GABAb_DynSyn_reg();
  _Glycine_DynSyn_reg();
  _HH2_reg();
  _HH2new_reg();
  _KDR_reg();
  _KDRI_reg();
  _NK1_DynSyn_reg();
  _NMDA_DynSyn_reg();
  _SS_reg();
  _borgka_reg();
  _iCaAN_reg();
  _iCaL_reg();
  _iKCa_reg();
  _iNaP_reg();
  _vecevent_reg();
  _vsource_reg();
}

#if defined(__cplusplus)
}
#endif
