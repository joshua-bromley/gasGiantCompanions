import os
import pandas as pd
import numpy as np
import radvel
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.special as ss
import os
import pylab as pl
from scipy import optimize
import corner
import pickle
#import mpfit
 
 
def missed_pl():
 
    sample ='all_metalrich'
 
    plname = np.array(['Teegarden','LSPMJ21160234','HD285968','GJ3293','GJ251','GJ229A','Wolf1061','Wolf1069','YZCet','TOI1266','TOI1452','TOI1468','TOI1695','TOI2018','proxima_cen','Ross128','Ross508','TOI244','TOI663','L9859','L36338','LHS1140','LHS1815','LTT1445A','K2-25','K2-415','Kapteyn','K2-18','Kepler138','HIP22627','HIP54373','HIP83043','HNLib','K2-3','GJ4276','GJ9689','GJ9827','HD180617','HD238090','GJ3634','GJ3779','GJ3929','GJ3988','GJ3998','GJ3082','GJ3138','GJ3323','GJ3341','GJ3470','GJ1002','GJ1132','GJ1151','GJ1214','GJ1265','GJ724','GJ740','GJ806','GJ876','GJ887','GJ674','GJ685','GJ686','GJ687','GJ720A','GJ514','GJ536','GJ581','GJ625','GJ667C','GJ411','GJ422','GJ433','GJ480','GJ486','GJ338B','GJ357','GJ367','GJ378','GJ393','GJ15A','GJ27','GJ49','GJ163','GJ273','AUMic','BD-082823','CD_Cet','G264-012','GJ12'])
   
    mass_star = np.array([0.097,0.43,0.5,0.42,0.36,0.509,0.294,0.167,0.142,0.431,0.249,0.339,0.54,0.57,0.1221,0.168,0.1774,0.428,0.514,0.273,0.21,0.179,0.502,0.257,0.2634,0.1635,0.281,0.359,0.535,0.393,0.57,0.52,0.291,0.549,0.406,0.59,0.62,0.484,0.578,0.45,0.27,0.313,0.184,0.50,0.47,0.68,0.164,0.47,0.539,0.12,0.181,0.164,0.182,0.178,0.527,0.58,0.413,0.370,0.489,0.35,0.55,0.426,0.40,0.57,0.51,0.52,0.311,0.3,0.33,0.39,0.35,0.48,0.45,0.323,0.64,0.342,0.455,0.56,0.426,0.38,0.53,0.515,0.38,0.29,0.51,0.50,0.161,0.297,0.241])
   
    metal_star = np.array([-0.11,-0.05,0.15,0.02,-0.03,0.15,-0.09,0.07,-0.18,-0.20,-0.07,-0.04,0.,-0.58,-0.04,-0.02,-0.20,-0.39,0.07,-0.46,-0.16,-0.15,-0.12,-0.34,0.15,-0.13,-0.89,0.12,-0.18,0.30,0.03,-0.15,-0.18,-0.16,0.12,0.05,-0.04,-0.03,-0.05,0.,-0.02,-0.12,-0.16,-0.22,0.,-0.27,-0.09,0.20,-0.25,-0.12,-0.12,0.24,-0.04,-0.02,0.08,-0.28,0.21,-0.06,-0.28,-0.15,-0.23,0.06,-0.14,-0.14,-0.08,-0.09,-0.38,-0.55,-0.36,0.18,-0.22,-0.06,0.07,0.,-0.12,-0.01,0.06,-0.09,-0.34,-0.09,0.13,-0.02,0.09,0.12,0.,0.13,0.10,-0.30])
 
    sel = np.where(plname == 'GJ3138')
    plname = np.delete(plname,sel)
    mass_star = np.delete(mass_star,sel)
    metal_star = np.delete(metal_star,sel)
    sel2 = np.where(plname == 'GJ338B')
    plname = np.delete(plname,sel2)
    mass_star = np.delete(mass_star,sel2)
    metal_star = np.delete(metal_star,sel2)
    sel3 = np.where(plname == 'BD-082823')
    plname = np.delete(plname,sel3)
    mass_star = np.delete(mass_star,sel3)
    metal_star = np.delete(metal_star,sel3)
    sel4 = np.where(plname == 'GJ9827')
    plname = np.delete(plname,sel4)
    mass_star = np.delete(mass_star,sel4)
    metal_star = np.delete(metal_star,sel4)
   
 
 
 
    sel_hm = np.where(metal_star > 0)
    sel_lm = np.where(metal_star <= 0)
 
    plname_metalrich = plname[sel_hm]
    plname_metalpoor = plname[sel_lm]
    mass_metalrich = mass_star[sel_hm]
    mass_metalpoor = mass_star[sel_lm]
    print('num metal rich')
    print(len(plname_metalrich))
    print('num metal poor')
    print(len(plname_metalpoor))
 
    #break the low and high metallicity samples into low and high stellar masses (0.3 M_sun is the cutoff)
    sel_hmhm = np.where(mass_metalrich > 0.3)
    sel_hmlm = np.where(mass_metalrich <=0.3)
    sel_lmhm = np.where(mass_metalpoor > 0.3)
    sel_lmlm = np.where(mass_metalpoor <= 0.3)
    plname_hmethmass = plname_metalrich[sel_hmhm]
    plname_hmetlmass = plname_metalrich[sel_hmlm]
    plname_lmethmass = plname_metalpoor[sel_lmhm]
    plname_lmetlmass = plname_metalpoor[sel_lmlm]
    print('num metal rich and high mass')
    print(len(plname_hmethmass))
    print('num metal rich and low mass')
    print(len(plname_hmetlmass))
    print('num metal poor and high mass')
    print(len(plname_lmethmass))
    print('num metal poor andlow mass')
    print(len(plname_lmetlmass))
   
 
    #break metallicity up into three bins
 
    sel_low = np.where(metal_star < -0.1)
    sel_med = np.where((metal_star>=-0.1)&(metal_star<=0.1))
    sel_high = np.where(metal_star > 0.1)
 
    plname_low = plname[sel_low]
    plname_med = plname[sel_med]
    plname_high = plname[sel_high]
    print('num low metal')
    print(len(plname_low))
    print('num mid metal')
    print(len(plname_med))
    print('num high metal')
    print(len(plname_high))
 
 
    if sample == 'all_metalrich':
        missed_pl = np.zeros(len(plname_lmethmass))
        for x in range(len(plname_lmethmass)):
            sys = plname_lmethmass[x]
            comp_prob = pickle.load(open('complete_maps_mstar_paper/'+sys+'_complete_map.p','rb'))
            probs = comp_prob['prob']
 
            sem_a_min = np.log10(0.3)
            sem_a_max = np.log10(30.)
            mass_min = np.log10(0.3)
            mass_max = np.log10(30.)
   
            exp_a = np.linspace(sem_a_min,sem_a_max,50)#np.arange(sem_a_min,sem_a_max,step_a)
            exp_m = np.linspace(mass_min,mass_max,50)#np.arange(mass_min,mass_max,step_m)
            val_a = 10.**exp_a
            sem_a = val_a
            val_m = 10.**exp_m
            mass = val_m
 
            low_a= find_nearest(sem_a,0.1)#1.)#(sem_a,0.1)#
            high_a = find_nearest(sem_a,5.)#10.)
            low_m = find_nearest(mass,0.5)#(mass,0.3)#
            high_m = find_nearest(mass,20.)
            sel_lowa = np.where(sem_a == low_a)
            #print(sel_lowa)
            sel_higha = np.where(sem_a == high_a)
            #print(sel_higha)
            sel_lowm = np.where(mass == low_m)
            #print(sel_lowm)
            sel_highm = np.where(mass == high_m)
            #print(sel_highm)
 
            #total could have been if 100% detected
            tot = (37.*40.)#(30.*40.)#
            sum_prob = np.sum(probs[0:37,5:45])#[0:30,5:45])#
            missed_pl[x] = (tot-sum_prob)/tot
            #print(missed_pl[x])
 
        print('total miss')
        print(missed_pl)
        total_miss = np.sum(missed_pl)
        print(total_miss)
      
    
 
 
def find_nearest(array,value):
    array = np.asarray(array)
    idx = (np.abs(array-value)).argmin()
    return array[idx]