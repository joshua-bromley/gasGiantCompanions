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
#import mpfit
 
 
def complete_plots():
    import pickle
    average = 'yes'#'minmax'#'minmax_total'#'nogg_median'#'no' #'gg'#'nogg'
    #average_trend = 'no'
 
    mcomp = np.array([0.752,2.108,0.698])
    semacomp = np.array([2.52,0.22,0.14])
 
    if average == 'yes':
        #plname = np.array(['Teegarden','LSPMJ21160234','HD285968','GJ3293','GJ251','GJ229A','Wolf1061','Wolf1069','YZCet','TOI1266','TOI1452','TOI1468','TOI1695','TOI2018','proxima_cen','Ross128','Ross508','TOI244','TOI663','L9859','L36338','LHS1140','LHS1815','LTT1445A','K2-25','K2-415','Kapteyn','K2-18','Kepler138','HIP22627','HIP54373','HIP83043','HNLib','K2-3','GJ4276','GJ9689','HD180617','HD238090','GJ3634','GJ3779','GJ3929','GJ3988','GJ3998','GJ3082','GJ3323','GJ3341','GJ3470','GJ1002','GJ1132','GJ1151','GJ1214','GJ1265','GJ724','GJ740','GJ806','GJ876','GJ887','GJ674','GJ685','GJ686','GJ687','GJ720A','GJ514','GJ536','GJ581','GJ625','GJ667C','GJ411','GJ422','GJ433','GJ480','GJ486','GJ357','GJ367','GJ378','GJ393','GJ15A','GJ27','GJ49','GJ163','GJ273','AUMic','CD_Cet','G264-012','GJ12'])
        plname = np.array(["Kepler407"])
        import pickle
        print(len(plname))
        prob = np.zeros((len(plname),49,49))
        for i in range(len(plname)):
            #print(i)
            comp_prob = pickle.load(open(f'./completenessMaps/{plname[i]}.p','rb'))
            prob[i,:,:] = comp_prob['prob']
            #print(prob[i,10,10])
 
        complete_ave = np.zeros((49,49))
        for j in range(49):
            for o in range(49):
                complete_ave[j,o] = np.sum(prob[:,j,o])
 
        print(complete_ave)
        prob = complete_ave/len(plname)
        sys='CLS_average'
       
    
    sem_a_min = np.log10(0.1)
    sem_a_max = np.log10(30.)
    mass_min = np.log10(0.3)
    mass_max = np.log10(30.)
   
    exp_a = np.linspace(sem_a_min,sem_a_max,(50*50))#np.arange(sem_a_min,sem_a_max,step_a)
    exp_m = np.linspace(mass_min,mass_max,(50*50))#np.arange(mass_min,mass_max,step_m)
    val_a = 10.**exp_a
    sem_a = val_a
    val_m = 10.**exp_m
    mass = val_m
 
    exp_a = np.linspace(sem_a_min,sem_a_max,(50))#np.arange(sem_a_min,sem_a_max,step_a)
    exp_m = np.linspace(mass_min,mass_max,(50))#np.arange(mass_min,mass_max,step_m)
    val_a = 10.**exp_a
    sem_as = val_a
    val_m = 10.**exp_m
    masss = val_m
 
    if average == 'minmax':
        cm = plt.cm.get_cmap('inferno')#('plasma')#('PuBu_r')
        massp,semap = np.meshgrid(mass[:-50],sem_a[:-50])
        masssp,semasp = np.meshgrid(masss[:-1],sem_as[:-1])
 
 
        plt.figure(figsize=(13,10))
        plt.rcParams['axes.linewidth']=2
        plt.xlabel('Semi-major axis [AU]',size=30)
        plt.ylabel('Mass [M'r'$\rm_{Jup}$'']',size=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.ylim([0.3,25])#375])
        plt.xlim([0.3,25])
        plt.yscale('log')
        plt.xscale('log')
        plt.contour(semasp,masssp,probmin,levels=[0.8],linewidths=3)
        plt.contour(semasp,masssp,probmax,levels=[0.8],linewidths=3)
        #plt.savefig('complete_maps_RVsys/'+sys+'_complete_maps.png')
        plt.show()
   
    # prob_s = np.zeros((len(mass)/2,len(sem_a)/2))
    # mass_n = np.zeros(len(mass)/2)
    # sem_a_n = np.zeros(len(sem_a)/2)
    #for i in xrange(len(sem_a)/2):
    #    for j in xrange(len(mass)/2):
    #        prob_s[i,j] = prob[i*2,j*2]
    #       mass_n[j] = mass[j*2]
    #       sem_a_n[i] = sem_a[i*2]
 
   
            
    #Mass_n,Sem_a_n, = np.meshgrid(mass_n,sem_a_n,)
    ones = np.ones((49,49))
    new_comp = np.zeros((49*50,49*50))
    for h in range(48):
        #print h
        for g in range(48):
            new_comp[(h*49+h):(h*49+h+49),(g*49+g):(g*49+g+49)] = ones*prob[h,g]#np.repeat(prob[h,g],500.,500.)
            new_comp_t = new_comp[(h*49+h):(h*49+h+49),(g*49+g):(g*49+g+49)]
           
 
    massp,semap = np.meshgrid(mass[:-50],sem_a[:-50])
    #print(massp.shape)
    #print(semap.shape)
    #print(new_comp.shape)
    masssp,semasp = np.meshgrid(masss[:-1],sem_as[:-1])
    #print(masssp.shape)
    #print(semasp.shape)
    #print(prob.shape)
 
    cm = plt.cm.get_cmap('inferno')#('plasma')#('PuBu_r')
 
    plt.figure(figsize=(15,10))
    plt.xlabel('Semi-Major Axis [AU]',size=30)
    plt.ylabel('Mass [M'r'$\rm_{Jup}$'']',size=30)
    #plt.rcParams['axes.linewidth']=2
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.ylim([0.3,25])#375])
    plt.xlim([0.1,23])
    plt.yscale('log')
    plt.xscale('log')
    plt.tick_params(width=2,length=10)
    #plt.contourf(semap,massp,new_comp,levels=[1.,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.])
    cont = plt.contourf(semap,massp,new_comp,levels=[0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],cmap=cm)#cmap = cm.PuBu_r)
    #plt.contourf(Sem_a_n,Mass_n,prob_s,levels=[1.,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.])#[0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.])#,colors=('blue','green','purple'),linewidths=(3,3,3))#,cmap = cm.PuBu_r)
    cbar = plt.colorbar(cont)
    cbar.ax.tick_params(labelsize=25)
    cbar.set_label('Detection Probability',rotation=270,fontsize=25,labelpad=25)
    #plt.tick_params(labelsize=15)
    # plt.set_label('Detection Probability',rotation=270,fontsize=20)
    #plt.plot(semacomp,mcomp,linestyle='none',marker='o',markersize=10,color='teal')
    #plt.savefig('complete_maps_rvsys_paper/'+sys+'_complete_maps.png')
    plt.savefig(f"./plots/{plname[i]}CompletenessMap.png")
    plt.show()

complete_plots()