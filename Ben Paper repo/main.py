# %% 
from neuron import h,gui
import matplotlib.pyplot as plt
import plotly
import plotly.graph_objects as go

# def __init__(self, L=45.75, nseg=1, diam=45.75, Ra=150.0, cm=1,
#                  pas_g=1.372e-5, pas_e=-58.0,
#                  HH2_gnabar=0.0126, HH2_gkbar=0.0286, HH2_vtraub=-32,
#                  KDR_gkbar= 1.6e-5,borgka_gkabar=2.29e-4,borgka_vhalfn=10,borgka_vhalfl=10,borgka_vhalfm=-49.32,
#                  borgka_vhalfk=-24.28,borgka_zetan=-0.1,borgka_zetal=0.1,borgka_gmn=0.5,borgka_gml=0.5,borgka_zetam=3.04,
#                  borgka_zetak=-2.54,borgka_a0n=0.012,borgka_a0l_slow=0.0018,borgka_a0l_fast=0.0096,borgka_fast_slow_ratio=0.62):
        

class PointNeuron:
    def __init__(self, L=45.75, nseg=1, diam=45.75, Ra=150.0, cm=1,
                 pas_g=1.372e-5, pas_e=-58.0,
                 HH2_gnabar=0.0126, HH2_gkbar=0.0286, HH2_vtraub=-32,
                 KDR_gkbar= 1.6e-5,borgka_gkabar=2.29e-4,borgka_vhalfn=10,borgka_vhalfl=10,borgka_vhalfm=-49.32,
                 borgka_vhalfk=-24.28,borgka_zetan=-0.1,borgka_zetal=0.1,borgka_gmn=0.5,borgka_gml=0.5,borgka_zetam=3.04,
                 borgka_zetak=-2.54,borgka_a0n=0.012,borgka_a0l_slow=0.0018,borgka_a0l_fast=0.0096,borgka_fast_slow_ratio=0.62):
        
        # Create neuron object according to specification and mechanisms above
        self.soma = h.Section(name='soma')
        self.soma.L = L
        self.soma.nseg = nseg
        self.soma.diam = diam
        self.soma.Ra = Ra
        self.soma.cm = cm

        # Insert mechanisms
        self.soma.insert('pas')
        self.soma.g_pas = pas_g
        self.soma.e_pas = pas_e

        self.soma.insert('HH2new')
        self.soma.gnabar_HH2new = HH2_gnabar
        self.soma.gkbar_HH2new = HH2_gkbar
        self.soma.vtraub_HH2new = HH2_vtraub
        

        self.soma.insert('KDR')
        self.soma.gkbar_KDR = KDR_gkbar

        self.soma.insert('borgka')
        self.soma.gkabar_borgka = borgka_gkabar

        self.soma.vhalfm_borgka = borgka_vhalfm
        self.soma.vhalfl_borgka = borgka_vhalfl
        self.soma.vhalfn_borgka = borgka_vhalfn
        self.soma.vhalfk_borgka = borgka_vhalfk

        self.soma.zetal_borgka = borgka_zetal
        self.soma.zetan_borgka=borgka_zetan
        self.soma.zetam_borgka=borgka_zetam
        self.soma.zetak_borgka=borgka_zetak

        self.soma.gmn_borgka=borgka_gmn
        self.soma.gml_borgka=borgka_gml

        self.soma.a0n_borgka=borgka_a0n
        self.soma.a0l_fast_borgka=borgka_a0l_fast
        self.soma.a0l_slow_borgka=borgka_a0l_slow
        self.soma.fsr_borgka=borgka_fast_slow_ratio
        
        #ek was set to match bens experiment
        self.set_ek(-88.85)

    def change_params(self, **params):
        for key, value in params.items():
            setattr(self, key, value)

    def set_ek(self,ek):
        for seg in self.soma:
            seg.ek=ek  
        
class PointNeuron_2:
    def __init__(self, L=20.0, nseg=1, diam=20.0, Ra=150.0, cm=6,
                 pas_g=5.15e-5, pas_e=-48,
                 HH2_gnabar=0.2, HH2_gkbar=0.4, HH2_vtraub=-32,
                 KDR_gkbar= 0e-5,borgka_gkabar=1.17e-3,borgka_vhalfn=-40,borgka_vhalfl=-40,borgka_vhalfm=-49.32,
                 borgka_vhalfk=-24.28,borgka_zetan=-1,borgka_zetal=1,borgka_gmn=0.5,borgka_gml=0.5,borgka_zetam=3.04,
                 borgka_zetak=-2.54,borgka_a0n=0.012,borgka_a0l=0.0024):
        
        # Create neuron object according to specification and mechanisms above
        self.soma = h.Section(name='soma')
        self.soma.L = L
        self.soma.nseg = nseg
        self.soma.diam = diam
        self.soma.Ra = Ra
        self.soma.cm = cm

        # Insert mechanisms
        self.soma.insert('pas')
        self.soma.g_pas = pas_g
        self.soma.e_pas = pas_e

        self.soma.insert('HH2new')
        self.soma.gnabar_HH2new = HH2_gnabar
        self.soma.gkbar_HH2new = HH2_gkbar
        self.soma.vtraub_HH2new = HH2_vtraub
        

        self.soma.insert('KDR')
        self.soma.gkbar_KDR = KDR_gkbar

        self.soma.insert('borgka')
        self.soma.gkabar_borgka = borgka_gkabar

        self.soma.vhalfm_borgka = borgka_vhalfm
        self.soma.vhalfl_borgka = borgka_vhalfl
        self.soma.vhalfn_borgka = borgka_vhalfn
        self.soma.vhalfk_borgka = borgka_vhalfk

        self.soma.zetal_borgka = borgka_zetal
        self.soma.zetan_borgka=borgka_zetan
        self.soma.zetam_borgka=borgka_zetam
        self.soma.zetak_borgka=borgka_zetak

        self.soma.gmn_borgka=borgka_gmn
        self.soma.gml_borgka=borgka_gml

        self.soma.a0n_borgka=borgka_a0n
        self.soma.a0l_borgka=borgka_a0l
        
        #ek was set to match bens experiment
        self.set_ek(-88.85)

    def change_params(self, **params):
        for key, value in params.items():
            setattr(self, key, value)

    def set_ek(self,ek):
        for seg in self.soma:
            seg.ek=ek  
        

class Soma:
    def __init__(self, L=20.0, nseg=1, diam=20.0, Ra=150.0, cm=1,
                 pas_g=3.7e-5, pas_e=-60.5,
                 HH2_gnabar=0.0, HH2_gkbar=0.001, HH2_vtraub=-55.0,
                 iCaL_pcabar=0.0001,
                 iCaAN_gbar=0.0,
                 iKCa_gbar=1e-5, iKCa_gk=0.0,
                 iNaP_gnabar=0.0001, iNaP_gamma=0.5, iNaP_vsh=-5.0, iNaP_vsm=-2.0, iNaP_vtraub=-55.0,
                 CaIntraCellDyn_cai_inf=5e-05, CaIntraCellDyn_cai_tau=1.0, CaIntraCellDyn_depth=0.1,
                 KDR_gkbar= 1e-5,borgka_gkabar=3.7e-4,borgka_vhalfn=-40,borgka_vhalfl=-40,borgka_vhalfm=-49.32,
                 borgka_vhalfk=-24.28,borgka_zetan=-1,borgka_zetal=0.1,borgka_gmn=0.5,borgka_gml=0.5,borgka_zetam=3.04,
                 borgka_zetak=-2.54,borgka_a0n=0.012,borgka_a0l=0.0024):
        # Create neuron object according to specification and mechanisms above
        self.soma = h.Section(name='soma')
        self.soma.L = L
        self.soma.nseg = nseg
        self.soma.diam = diam
        self.soma.Ra = Ra
        self.soma.cm = cm

        # Insert mechanisms
        self.soma.insert('pas')
        self.soma.g_pas = pas_g
        self.soma.e_pas = pas_e

        self.soma.insert('HH2')
        self.soma.gnabar_HH2 = HH2_gnabar
        self.soma.gkbar_HH2 = HH2_gkbar
        self.soma.vtraub_HH2 = HH2_vtraub

        # self.soma.insert('iCaL')
        # self.soma.pcabar_iCaL = iCaL_pcabar

        # self.soma.insert('iCaAN')
        # self.soma.gbar_iCaAN = iCaAN_gbar

        # self.soma.insert('iKCa')
        # self.soma.gbar_iKCa = iKCa_gbar
        # self.soma.gk_iKCa = iKCa_gk

        # self.soma.insert('iNaP')
        # self.soma.gnabar_iNaP = iNaP_gnabar
        # self.soma.gamma_iNaP = iNaP_gamma
        # self.soma.vsh_iNaP = iNaP_vsh
        # self.soma.vsm_iNaP = iNaP_vsm
        # self.soma.vtraub_iNaP = iNaP_vtraub

        # self.soma.insert('CaIntraCellDyn')
        # self.soma.cai_inf_CaIntraCellDyn = CaIntraCellDyn_cai_inf
        # self.soma.cai_tau_CaIntraCellDyn = CaIntraCellDyn_cai_tau
        # self.soma.depth_CaIntraCellDyn = CaIntraCellDyn_depth

        """kdr was removed from the model and replaced with borgka, the observed reting membrane potentials 
        were not consistent with the existance of a KDR current of a siginificant magnitude"""
        self.soma.insert('KDR')
        self.soma.gkbar_KDR = KDR_gkbar

        self.soma.insert('borgka')
        self.soma.gkabar_borgka = borgka_gkabar

        self.soma.vhalfm_borgka = borgka_vhalfm
        self.soma.vhalfl_borgka = borgka_vhalfl
        self.soma.vhalfn_borgka = borgka_vhalfn
        self.soma.vhalfk_borgka = borgka_vhalfk

        self.soma.zetal_borgka = borgka_zetal
        self.soma.zetan_borgka=borgka_zetan
        self.soma.zetam_borgka=borgka_zetam
        self.soma.zetak_borgka=borgka_zetak

        self.soma.gmn_borgka=borgka_gmn
        self.soma.gml_borgka=borgka_gml

        self.soma.a0n_borgka=borgka_a0n
        self.soma.a0l_borgka=borgka_a0l
        
        #ek was set to match bens experiment
        self.set_ek(-88.85)

    def change_params(self, **params):
        for key, value in params.items():
            setattr(self, key, value)

    def set_ek(self,ek):
        for seg in self.soma:
            seg.ek=ek  


class Dendrite:
    def __init__(self, L=350.0, nseg=5, diam=2.5, Ra=150.0, cm=1,
                 pas_g=3.7e-5, pas_e=-60.5,
                 HH2_gnabar=0.0, HH2_gkbar=0.01, HH2_vtraub=-35.0,
                 iCaL_pcabar=3e-05,
                 iCaAN_gbar=9.099999999999999e-05,
                 iKCa_gbar=5e-5, iKCa_gk=0.0,
                 CaIntraCellDyn_cai_inf=5e-05, CaIntraCellDyn_cai_tau=2.0, CaIntraCellDyn_depth=0.1,
                 KDR_gkbar= 5e-5):
        # Create neuron object according to specification and mechanisms above
        self.dend = h.Section(name='dendrite')
        self.dend.L = L
        self.dend.nseg = nseg
        self.dend.diam = diam
        self.dend.Ra = Ra
        self.dend.cm = cm

        # Insert mechanisms
        self.dend.insert('pas')
        self.dend.g_pas = pas_g
        self.dend.e_pas = pas_e

        self.dend.insert('HH2')
        self.dend.gnabar_HH2 = HH2_gnabar
        self.dend.gkbar_HH2 = HH2_gkbar
        self.dend.vtraub_HH2 = HH2_vtraub

        # self.dend.insert('iCaL')
        # self.dend.pcabar_iCaL = iCaL_pcabar

        # self.dend.insert('iCaAN')
        # self.dend.gbar_iCaAN = iCaAN_gbar

        # self.dend.insert('iKCa')
        # self.dend.gbar_iKCa = iKCa_gbar
        # self.dend.gk_iKCa = iKCa_gk

       
        # self.dend.insert('CaIntraCellDyn')
        # self.dend.cai_inf_CaIntraCellDyn = CaIntraCellDyn_cai_inf
        # self.dend.cai_tau_CaIntraCellDyn = CaIntraCellDyn_cai_tau
        # self.dend.depth_CaIntraCellDyn = CaIntraCellDyn_depth

        """kdr was removed from the model and replaced with borgka, the observed reting membrane potentials 
        were not consistent with the existance of a KDR current of a siginificant magnitude"""
        self.dend.insert('KDR')
        self.dend.gkbar_KDR = KDR_gkbar
        
        #ek was set to match bens experiment
        self.set_ek(-88.85)
        

    def change_params(self, **params):
        for key, value in params.items():
            setattr(self, key, value)

    def set_ek(self,ek):
        for seg in self.dend:
            seg.ek=ek       

        


class Hillock:
    def __init__(self, L=9.0, nseg=1, diam=1.5, Ra=150.0, cm=1,
                 pas_g=3.7e-5, pas_e=-60.5,
                 HH2_gnabar=1, HH2_gkbar=0.7, HH2_vtraub=-39.0):
        # Create neuron object according to specification and mechanisms above
        self.hil = h.Section(name='hillock')
        self.hil.L = L
        self.hil.nseg = nseg
        self.hil.diam = diam
        self.hil.Ra = Ra
        self.hil.cm = cm

        # Insert mechanisms
        self.hil.insert('pas')
        self.hil.g_pas = pas_g
        self.hil.e_pas = pas_e

        self.hil.insert('HH2')
        self.hil.gnabar_HH2 = HH2_gnabar
        self.hil.gkbar_HH2 = HH2_gkbar
        self.hil.vtraub_HH2 = HH2_vtraub

        #ek was set to match bens experiment
        self.set_ek(-88.85)

    def change_params(self, **params):
        for key, value in params.items():
            setattr(self, key, value)

    def set_ek(self,ek):
        for seg in self.hil:
            seg.ek=ek          


       

# Run the simulation
# h.dt = 0.025
# h.tstop = 1100
# h.celsius = 37
# h.v_init = -61
# v_vec = h.Vector()
# t_vec = h.Vector()



# soma=Soma()
# dend=Dendrite()
# hil=Hillock()


# soma.soma.connect(hil.hil(0))
# dend.dend.connect(soma.soma(1))




# # ps=h.PlotShape(False)
# # ps.variable("v")
# # ps.plot(plotly).show()

# # %% 



# # Add current stimulation
# stim = h.IClamp(soma.soma(0.5))
# stim.delay = 100
# stim.dur = 1000
# stim.amp = 0.00


# # # Record the membrane potential
# v_vec.record(soma.soma(0.5)._ref_v)
# t_vec.record(h._ref_t)

# # # # 
# h.run()

# # # # # %% Plot the results
# import matplotlib.pyplot as plt
# plt.plot(t_vec, v_vec)
# plt.show()

# %%

