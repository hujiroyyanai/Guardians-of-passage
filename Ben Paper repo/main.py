# %% 
from neuron import h,gui
import matplotlib.pyplot as plt
import plotly
import plotly.graph_objects as go


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
         
        


       

