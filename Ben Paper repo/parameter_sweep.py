from neuron import h,gui
import matplotlib as mpl
import matplotlib.pyplot as plt
from main import Soma,Dendrite, Hillock
from main import PointNeuron , PointNeuron_2
import pandas as pd
import numpy as np

#code for working with mods files
def transform_a0l_to_tau(a0l):
    return 1/(a0l*2*(3**0.7))
    
#code for working with mods files
def transform_tau_to_a0l(tau):
    return 1/(tau*2*(3**0.7))



def check_validity(simulated_value,valid_range):
    if simulated_value < valid_range[0] or simulated_value > valid_range[1]:
        return False
    else:
        return True



def validate_firing_charactaristics(cell,amp,latency_valid_range,frequency_valid_range):
    h.dt = 0.025
    h.tstop = 1100
    h.celsius = 37
    h.v_init = -62
    v_vec = h.Vector()
    t_vec = h.Vector()
    spike_times = h.Vector()


    # Add current stimulation
    stim = h.IClamp(cell.soma(0.5))
    stim.delay = 100
    stim.dur = 500
    stim.amp = amp

    # Record the membrane potential
    v_vec.record(cell.soma(0.5)._ref_v)
    t_vec.record(h._ref_t)


    #Record spike times
    nc=h.NetCon(cell.soma(0.5)._ref_v,None,sec=cell.soma)
    nc.threshold = 0
    nc.record(spike_times)

    # 
    h.run()

    import matplotlib.pyplot as plt
    plt.plot(t_vec, v_vec)
    plt.show()
    
    print(2*len(spike_times.x))

    if len(spike_times.x) == 0:
        return False
    
    if not (spike_times.x[0] > latency_valid_range[0] and spike_times.x[0] < latency_valid_range[1]):
        return False
    
    #current was injected for 500 ms, so the frequency is the number of spikes divided by 0.5
    frequency_of_fire=len(spike_times.x)/0.5 

    if not (frequency_of_fire > frequency_valid_range[0] and frequency_of_fire < frequency_valid_range[1]):
        return False
    
    return True

def validate_resting_potential(cell,valid_range,inflamation=False):
    if inflamation==True:
        cell.soma.vhalfm_borgka=-45.2079

    h.dt = 0.25
    h.tstop = 1100
    h.celsius = 37
    h.v_init = -62
    v_vec = h.Vector()
    t_vec = h.Vector()

    # Record the membrane potential
    v_vec.record(cell.soma(0.5)._ref_v)
    t_vec.record(h._ref_t)

    # # 
    h.run()
    resting_potential = v_vec.x[-1]
    # plt.plot(t_vec, v_vec)
    # plt.show()

    return (check_validity(resting_potential,valid_range),resting_potential)

def validate_input_resistance(cell,resting_potential,valid_range, inflamation=False):
    if inflamation==True:
        cell.soma.vhalfm_borgka=-45.2079

    h.dt = 0.25
    h.tstop = 1100
    h.celsius = 37
    h.v_init = -62
    v_vec = h.Vector()
    t_vec = h.Vector()

    # Record the membrane potential
    v_vec.record(cell.soma(0.5)._ref_v)
    t_vec.record(h._ref_t)

    # # 
    stim = h.IClamp(cell.soma(0.5))
    stim.delay = 100
    stim.dur = 1000
    stim.amp = -0.001

    h.run()

    delta_v = v_vec.x[-1] - resting_potential
    #R = delta_v/I in Momhs
    resistance=(1e-3*delta_v/(1e-9*stim.amp))*1e-6
    # plt.plot(t_vec, v_vec)
    # plt.show()

    return (check_validity(resistance,valid_range),resistance)

def calc_curves(cell,amp_range):
    """clac IF and letancy curves relative to current amplitude"""
    fig, axs = plt.subplots(2, 2)
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    f=[]
    l=[]
    for amp in amp_range:
        h.dt = 0.025
        h.tstop = 1100
        h.celsius = 37
        h.v_init = -62
        v_vec = h.Vector()
        t_vec = h.Vector()
        spike_times = h.Vector()

        stim = h.IClamp(cell.soma(0.5))
        stim.delay = 100
        stim.dur = 500
        stim.amp = amp

        # Record the membrane potential
        v_vec.record(cell.soma(0.5)._ref_v)
        t_vec.record(h._ref_t)

        #Record spike times
        nc=h.NetCon(cell.soma(0.5)._ref_v,None,sec=cell.soma)
        nc.threshold = 0
        nc.record(spike_times)

        # 
        h.run()

        freq=2*len(spike_times.x)
        f.append(freq)

        if freq>0:
            if spike_times.x[0]<140:
                l.append(spike_times.x[1]-100)
            else:
                l.append(spike_times.x[0]-100)    
    
        else:
            l.append(-100)    

        if amp==0.08 :
            axs[0,0].plot(t_vec, v_vec)
            axs[0,0].set_title('Amp=0.08')
            
        elif amp==0.1:
            axs[1,0].plot(t_vec, v_vec)  
            axs[1,0].set_title('Amp=0.1') 

    axs[0,1].scatter(amp_range,f)
    axs[0,1].set_title('Firing frequency')

    axs[1,1].scatter(amp_range,l)
    axs[1,1].set_title('Latency')
    plt.show()


def calc_curves_with_inflamation(cell,amp_range):
    """clac IF and letancy curves relative to current amplitude"""
    data_dict={'frequencies':[None,None],'latancies':[None,None],'inst_frequencies':[None,None]}

    fig, axs = plt.subplots(3, 2)
    plt.subplots_adjust(hspace=0.5, wspace=0.5)

    for inflamation in range(2):
        if inflamation==1:
            cell.soma.vhalfm_borgka=-45.2079

        f=[]
        l=[]
        inst_f=[]
        for i,amp in enumerate(amp_range):
            h.dt = 0.025
            h.tstop = 1010
            h.celsius = 37
            h.v_init = -70
            v_vec = h.Vector()
            t_vec = h.Vector()
            spike_times = h.Vector()

            voltage_clamp=h.SEClamp(cell.soma(0.5))
            voltage_clamp.dur1=100
            voltage_clamp.amp1=-70

            stim = h.IClamp(cell.soma(0.5))
            stim.delay = 100
            stim.dur = 500
            stim.amp = amp

            # Record the membrane potential
            v_vec.record(cell.soma(0.5)._ref_v)
            t_vec.record(h._ref_t)

            #Record spike times
            nc=h.NetCon(cell.soma(0.5)._ref_v,None,sec=cell.soma)
            nc.threshold = 0
            nc.record(spike_times)

            # 
            h.run()

            freq=2*len(spike_times.x)
            f.append(freq)
            
            #append latency to first spike
            if freq>0:
                if spike_times.x[0]<140:
                    l.append(spike_times.x[1]-100)
                else:
                    l.append(spike_times.x[0]-100)  

            else:
                l.append(float('nan'))        

            #compute instantanous freq
            num_spikes=len(spike_times.x)    
            if num_spikes>1:
                isi=[]
                for indx,spike in enumerate(spike_times):
                    if indx<num_spikes-1:
                        isi.append((spike_times.x[indx+1]-spike_times.x[indx])/1000)     

                isi_inverse=[1/el for el in isi]
                mean=sum(isi_inverse)/len(isi_inverse)
                inst_f.append(mean)
            else:
                inst_f.append(float('nan'))                
                

            if 0.08<=amp and amp<=0.12 :
                axs[i-1,inflamation].plot(t_vec, v_vec)
                axs[i-1,inflamation].set_title(f'Amp={amp}')

                axs[i-1,inflamation].set_xlim(0,1000)
                axs[i-1,inflamation].set_ylim(-80,40)
                axs[i-1,inflamation].set_xticks([0,200,400,600,800,1000])
                axs[i-1,inflamation].set_yticks([-80,-60,-40,-20,0,20,40])
    
        print(f"inst freq are {inst_f}")
        data_dict['frequencies'][inflamation]=f
        data_dict['latancies'][inflamation]=l
        data_dict['inst_frequencies'][inflamation]=inst_f
    plt.show()
    return data_dict


def calc_vhalf_param_sweep(cell,amp_range,vhalf_params):
    latancies=pd.DataFrame(np.nan,columns=amp_range,index=vhalf_params)
    frequencies=pd.DataFrame(np.nan,columns=amp_range,index=vhalf_params)

    for amp in amp_range:
        for vhalf in vhalf_params:
            cell.soma.vhalfm_borgka=vhalf
            h.dt = 0.025
            h.tstop = 1100
            h.celsius = 37
            h.v_init = -70
            v_vec = h.Vector()
            t_vec = h.Vector()
            spike_times = h.Vector()

            voltage_clamp=h.SEClamp(cell.soma(0.5))
            voltage_clamp.dur1=100
            voltage_clamp.amp1=-70

            stim = h.IClamp(cell.soma(0.5))
            stim.delay = 100
            stim.dur = 500
            stim.amp = amp

            # Record the membrane potential
            v_vec.record(cell.soma(0.5)._ref_v)
            t_vec.record(h._ref_t)

            #Record spike times
            nc=h.NetCon(cell.soma(0.5)._ref_v,None,sec=cell.soma)
            nc.threshold = 0
            nc.record(spike_times)

            # 
            h.run()

            freq=2*len(spike_times.x)

            if len(spike_times.x)>0:
                latancy=(spike_times.x[0]-100)
            else:   
                latancy=500 

            latancies.at[vhalf,amp]=latancy
            frequencies.at[vhalf,amp]=freq
    print(frequencies)

    fig,axs=plt.subplots(2)
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    axs[0].pcolor(frequencies,cmap='Reds_r')
    axs[0].set_title('Frequencies')
    axs[0].set_yticklabels(vhalf_params)
    axs[0].set_xticklabels(amp_range)
    fig.colorbar(axs[0].pcolor(frequencies,cmap='Reds_r'),label='Frequency (Hz)')
    axs[1].pcolor(latancies,cmap='Reds_r')
    axs[1].set_title('Latancies')
    axs[1].set_yticklabels(vhalf_params)
    axs[1].set_xticklabels(amp_range)
    fig.colorbar(axs[1].pcolor(latancies,cmap='Reds_r'),label='Latancy (ms)')
    plt.show()
    

            

    
 


def plot_IF_and_latancy(data_dict,amp_range):
    fig, axs = plt.subplots(3)
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    axs[0].scatter(x=amp_range,y=data_dict['frequencies'][0],color='blue',alpha=0.5)
    axs[0].plot(amp_range,data_dict['frequencies'][0],color='blue')
    axs[0].scatter(x=amp_range,y=data_dict['frequencies'][1],color='red',alpha=0.5)
    axs[0].plot(amp_range,data_dict['frequencies'][1],color='red')
    axs[0].set_title('FI')

    axs[1].scatter(x=amp_range,y=data_dict['inst_frequencies'][0],color='blue',alpha=0.5)
    axs[1].plot(amp_range,data_dict['inst_frequencies'][0],color='blue')
    axs[1].scatter(x=amp_range,y=data_dict['inst_frequencies'][1],color='red',alpha=0.5)
    axs[1].plot(amp_range,data_dict['inst_frequencies'][1],color='red')
    axs[1].set_title('inst_FI')

    axs[2].scatter(x=amp_range,y=data_dict['latancies'][0],color='blue')
    axs[2].plot(amp_range,data_dict['latancies'][0],color='blue')
    axs[2].scatter(x=amp_range,y=data_dict['latancies'][1],color='red')
    axs[2].plot(amp_range,data_dict['latancies'][1],color='red')
    axs[2].set_title('latancy')
    plt.show()

#these lines of code are for running the simulations presented in the paper  
toy_cell=PointNeuron()


_,resting_potential=validate_resting_potential(toy_cell,(-63,-60),inflamation=False)
print(resting_potential)

res=validate_input_resistance(toy_cell,resting_potential,(615-32,615+32),inflamation=False)
print(res)


data_dict=calc_curves_with_inflamation(toy_cell,[0.06,0.08,0.1,0.12,0.14])
plot_IF_and_latancy(data_dict,[0.06,0.08,0.1,0.12,0.14])



calc_vhalf_param_sweep(toy_cell,[0.06,0.08,0.1,0.12,0.14],np.linspace(-49.32,-45.207,5))







    


    
      

