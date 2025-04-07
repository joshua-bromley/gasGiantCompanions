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
 
 
def complete_plots(plname, masses, smaxes):
    import pickle
    average = 'yes'#'minmax'#'minmax_total'#'nogg_median'#'no' #'gg'#'nogg'
    #average_trend = 'no'
 
    #mcomp = np.array([1.981,4.132,9.45,16.35,1.40642,3.18,7.79,6.97,3.71,1.31,0.59,15,1.99,11.1,12.7,10.7,3.4,14.28,17])
    #mcomp = np.array([1.37,4.89,0.57,5.61,0.784,1.5,0.549,0.82,1.26,0.89,6.493,3.86,2.07659,1.69,1.249,10.901,1.33,1.18,0.66,0.622,4.8,12.8,0.93,4,2.51,3.12])
    #semacomp = np.array([0.827774,2.51329,5.28,3.73,1.403666,1.07,5.5,0.73,2.39,0.54,1.94,13.2,4.89,6.2,2.504,8.968117,5.6,1.25,2.9])
    #semacomp = np.array([0.4756,2,0.1652,2.16,0.316,8.00,0.902,2.07,1.074,1.955,5.362,4.92,3.2,1.14,1.039,9.811,0.689,2.09,5,2.16,2.885,18.2,2.44,1.7520,4.11,5.9570])
    mcomp = masses
    semacomp = smaxes

    if average == 'yes':
        #plname = np.array(['Teegarden','LSPMJ21160234','HD285968','GJ3293','GJ251','GJ229A','Wolf1061','Wolf1069','YZCet','TOI1266','TOI1452','TOI1468','TOI1695','TOI2018','proxima_cen','Ross128','Ross508','TOI244','TOI663','L9859','L36338','LHS1140','LHS1815','LTT1445A','K2-25','K2-415','Kapteyn','K2-18','Kepler138','HIP22627','HIP54373','HIP83043','HNLib','K2-3','GJ4276','GJ9689','HD180617','HD238090','GJ3634','GJ3779','GJ3929','GJ3988','GJ3998','GJ3082','GJ3323','GJ3341','GJ3470','GJ1002','GJ1132','GJ1151','GJ1214','GJ1265','GJ724','GJ740','GJ806','GJ876','GJ887','GJ674','GJ685','GJ686','GJ687','GJ720A','GJ514','GJ536','GJ581','GJ625','GJ667C','GJ411','GJ422','GJ433','GJ480','GJ486','GJ357','GJ367','GJ378','GJ393','GJ15A','GJ27','GJ49','GJ163','GJ273','AUMic','CD_Cet','G264-012','GJ12'])
        #plname = np.array(["upsAnd", "WASP8", "WASP53", "WASP47", "WASP41", "Pr0211", "Kepler424", "KELT6", "HIP14810", "HD68988", "HD187123", "HD118203", "HATS59", "HATP2", "HATP11", "HATP13", "CoRoT20"])
        #plname = np.array(['14Her', '24Sex', '47UMa', '55Cnc', '7CMa', '75Cet', 'BD202457', 'GJ317', 'GJ676' ,'GJ849', 'GJ876', 'HD102329', 'HD108874','HD111232' ,'HD113337', 'HD114783' ,'HD11506' ,'HD116029', 'HD125612','HD12661' ,'HD128311' ,'HD134987' ,'HD13908' ,'HD141399', 'HD142' ,'HD147873', 'HD148164', 'HD154857', 'HD155358', 'HD156279','HD159243', 'HD159868', 'HD1605', 'HD160691', 'HD168443' ,'HD169830','HD181433', 'HD183263', 'HD202696' ,'HD204313' ,'HD207832','HD27894' ,'HD30177' ,'HD33142', 'HD33844', 'HD37124', 'HD38529', 'HD4203', 'HD4732', 'HD47366', 'HD50499', 'HD5319', 'HD67087', 'HD72659', 'HD74156', 'HD75784','HD75898', 'HD92788', 'HD9446', 'HD99706', 'HIP14810','HIP67851', 'HIP8541', 'Kepler432', 'Kepler454', 'Kepler48','Kepler539' ,'Kepler56' ,'TIC279401253', 'TYC14226141','gamLib', 'upsAnd'])
        #plname = np.array(["55Cnc", "GJ328", "HATP44", "HD113538", "HD11506", "HD11964", "HD141399", "HD142", "HD177830", "HD191939", "HD24040", "HD27894", "HD33142", "HD83443", "HIP65407", "Kepler56","TIC139270665", "XO2S"])
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
    plt.plot(semacomp,mcomp,linestyle='none',marker='o',markersize=10,color='teal')
    #plt.savefig('complete_maps_rvsys_paper/'+sys+'_complete_maps.png')
    plt.savefig(f"./plots/hjCompanionsCompletenessMap.png")
    plt.show()


planets = pd.read_csv("./data/gasGiantDataComplete.csv")
planets = planets.loc[(pd.isna(planets["pl_bmassj"]) == False) & (pd.isna(planets["pl_orbsmax"]) == False) & (planets["pl_orbeccen"] > 0)]
coldJupiters = planets.loc[(planets["pl_type"] == "WJ") | (planets["pl_type"] == "CJ")]
seCompanions = planets.loc[planets["companion_type"] %2 == 0]
ssCompanions = planets.loc[(planets["companion_type"] % 3 == 0) | (planets["companion_type"] % 5 == 0)]
hjCompanions = planets.loc[planets["companion_type"] % 7 == 0]
cjCompanions = planets.loc[(planets["companion_type"] % 11 == 0) | (planets["companion_type"] % 13 == 0)]


hostnames = np.array(hjCompanions["hostname"])
masses = np.array(hjCompanions["pl_bmassj"])
smaxes = np.array(hjCompanions["pl_orbsmax"])
complete_plots(hostnames, masses, smaxes)