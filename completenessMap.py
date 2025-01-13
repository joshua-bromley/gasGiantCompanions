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
import radvel.likelihood
import corner
from PyAstronomy.pyTiming import pyPeriod
 
#import mpfit
 
 
def complete_maps():
 
    sys = np.array(['Teegarden'])#(['LSPMJ21160234'])#(['HD285968'])#(['GJ3293'])#(['GJ251'])#(['GJ229A'])#(['Wolf1061','Wolf1069','YZCet'])#(['TOI1266','TOI1452','TOI1468','TOI1695','TOI2018'])#(['proxima_cen','Ross128','Ross508','TOI244','TOI663'])#(['L9859','L36338','LHS1140','LHS1815','LTT1445A'])#(['K2-25','K2-415','Kapteyn','K2-18','Kepler138'])#(['HIP22627','HIP54373','HIP83043','HNLib','K2-3'])#(['GJ4276','GJ9689','GJ9827','HD180617','HD238090'])#(['GJ3634','GJ3779','GJ3929','GJ3988','GJ3998'])#(['GJ3082','GJ3138','GJ3323','GJ3341','GJ3470'])#(['GJ1002','GJ1132','GJ1151','GJ1214','GJ1265'])#(['GJ724','GJ740','GJ806','GJ876','GJ887'])#(['GJ674','GJ685','GJ686','GJ687','GJ720A'])#(['GJ514','GJ536','GJ581','GJ625','GJ667C'])#(['GJ411','GJ422','GJ433','GJ480','GJ486'])#(['GJ338B','GJ357','GJ367','GJ378','GJ393'])#(['GJ15A','GJ27','GJ49','GJ163','GJ273'])#(['AUMic','BD-082823','CD_Cet','G264-012','GJ12'])
   
    mass_star = np.array([0.097])#([0.43])#([0.5])#([0.42])#([0.36])#([0.509])#([0.294,0.167,0.142])#([0.431,0.249,0.339,0.54,0.57])#([0.1221,0.168,0.1774,0.428,0.514])#([0.273,0.21,0.179,0.502,0.257])#([0.2634,0.1635,0.281,0.359,0.535])#([0.393,0.57,0.52,0.291,0.549])#([0.406,0.59,0.62,0.484,0.578])#([0.45,0.27,0.313,0.184,0.50])#([0.47,0.68,0.164,0.47,0.539])#([0.12,0.181,0.164,0.182,0.178])#([0.527,0.58,0.413,0.370,0.489])#([0.35,0.55,0.426,0.40,0.57])#([0.51,0.52,0.311,0.3,0.33])#([0.39,0.35,0.48,0.45,0.323])#([0.64,0.342,0.455,0.56,0.426])#([0.38,0.53,0.515,0.38,0.29])#([0.51,0.50,0.161,0.297,0.241])
   
    for x in range(len(sys)):
        print(sys[x])
        plname = sys[x]
        n_inject = 50
 
        sem_a_min = np.log10(0.3)
        sem_a_max = np.log10(30.)
        mass_min = np.log10(0.3)
        mass_max = np.log10(30.)
 
        exp_a = np.linspace(sem_a_min,sem_a_max,50)
        exp_m = np.linspace(mass_min,mass_max,50)
        val_a = 10.**exp_a
        sem_a = val_a
        val_m = 10.**exp_m
        mass = val_m
 
       
        #For BD-08 2823, GJ 367, GJ 674, GJ 1132, GJ 3634, GJ 9827, LHS1140, Ross 128
        if plname == 'BD-082823' or plname == 'GJ367' or plname == 'GJ674' or plname == 'GJ1132' or plname == 'GJ3634' or plname == 'GJ9827' or plname == 'LHS1140' or plname == 'Ross128':
            time,rv,err = np.loadtxt('RV_datasets_mstar/RVs_standard/'+plname+'_rvs_std.csv',delimiter=',',unpack=True)
            time = time - 2450000
            err=err*1000.
        #------------------------------------
 
        #For HD15337, HD 73583, K2-180 WASP-47, K2-100, pi Men
        #if plname == 'HD15337' or plname == 'HD73583' or plname == 'K2180' or plname == 'K2216' or plname == 'K2100' or plname == 'piMen':
         #   time,rv,err = np.loadtxt('RV_datasets_transitsys/RVs_standard/'+plname+'_rvs_std.csv',delimiter=',',unpack=True)
          #  err=err*1000.
        #------------------------------------
 
        #For GJ 163, GJ 273, GJ 378, GJ 3138, GJ 3323, GJ 3341, TOI-663, Wolf 1061, GJ 3293
        if plname == 'GJ163' or plname == 'GJ273' or plname == 'GJ378' or plname == 'GJ3138' or plname == 'GJ3323' or plname == 'GJ3341' or plname == 'TOI663' or plname == 'Wolf1061' or plname == 'GJ3293':
            time,rv,err = np.loadtxt('RV_datasets_mstar/RVs_standard/'+plname+'_rvs_std.csv',delimiter=',',unpack=True)
            time = time - 50000
            err=err*1000.
        #------------------------------------
       
 
        #For AU Mic, CD Cet, G 264-012, GJ 12, GJ 15A, GJ 357, GJ 393, GJ 411, GJ 480, GJ 486, GJ 581, GJ 667C. GJ 686, GJ 687, GJ 724, GJ 740, GJ 806, GJ 876, GJ 887, GJ 1151, GJ 1265, GJ 3470, GJ 3779, GJ 3929, GJ 3988, GJ 4276, HD 180617, HD 238090, HIP 22627, HIP 54373, HIP 83043, HN Lib, K2-25, K2-415, Kapteyn, Kepler-138, LHS 1815, LTT 1445A, Proxima Cen, TOI-1452, TOI-1468, TOI-1695, TOI-2018, Wolf 1069, YZ Cet, GJ 251, HD 285968, LSPMJ21160234, Teegarden
        if plname == 'AUMic' or plname == 'CD_Cet' or plname == 'G264-012' or plname == 'GJ12' or plname == 'GJ15A' or plname == 'GJ357' or plname == 'GJ393' or plname == 'GJ411' or plname == 'GJ480' or plname == 'GJ486' or plname == 'GJ581' or plname == 'GJ667C' or plname == 'GJ686' or plname == 'GJ687' or plname == 'GJ724' or plname == 'GJ740' or plname == 'GJ806' or plname == 'GJ876' or plname == 'GJ887' or plname == 'GJ1151' or plname == 'GJ1265' or plname == 'GJ3470' or plname == 'GJ3779' or plname == 'GJ3929' or plname == 'GJ3988' or plname == 'GJ4276' or plname == 'HD180617' or plname == 'HD238090' or plname == 'HIP22627' or plname == 'HIP54373' or plname == 'HIP83043' or plname == 'HNLib' or plname == 'K2-25' or plname == 'K2-415' or plname == 'Kapteyn' or plname == 'Kepler138' or plname == 'LHS1815' or plname == 'LTT1445A' or plname == 'proxima_cen' or plname == 'TOI1452' or plname == 'TOI1468' or plname == 'TOI1695' or plname == 'TOI2018' or plname == 'Wolf1069' or plname == 'YZCet' or plname == 'GJ251' or plname == 'HD285968' or plname == 'LSPMJ21160234' or plname == 'Teegarden':
            time,rv,err = np.loadtxt('RV_datasets_mstar/RVs_standard/'+plname+'_rvs_std.csv',delimiter=',',unpack=True)
            time = time - 2450000#2450000
        #------------------------------------
 
        #For GJ 27, GJ 685, GJ 720A, GJ 9689, L 98-59, TOI-244
        if plname == 'GJ27' or plname == 'GJ685' or plname == 'GJ720A' or plname == 'GJ9689' or plname == 'L9859' or plname == 'TOI244':
            time,rv,err = np.loadtxt('RV_datasets_mstar/RVs_standard/'+plname+'_rvs_std.csv',delimiter=',',unpack=True)
            time = time-50000
 
        #For GJ 49, GJ 338B, GJ 422, GJ 433, GJ 514, GJ 536, GJ 625, GJ 1002, GJ 1214, GJ 3082, GJ 3998, K2-3, K2-18, L363-38, Ross 508, TOI-1266, GJ 229A
        if plname == 'GJ49' or plname == 'GJ338B' or plname == 'GJ422' or plname == 'GJ433' or plname == 'GJ514' or plname == 'GJ536' or plname == 'GJ625' or plname == 'GJ1002' or plname == 'GJ1214' or plname == 'GJ3082' or plname == 'GJ3998' or plname == 'K2-3' or plname == 'K2-18' or plname == 'L36338' or plname == 'Ross508' or plname == 'TOI1266' or plname == 'GJ229A':
            time,rv,err = np.loadtxt('RV_datasets_mstar/RVs_standard/'+plname+'_rvs_std.csv',delimiter=',',unpack=True)
 
       
 
        #------------------------------------
        jitter = 0.
        mstar = mass_star[x]
        print('mstar')
        print(mstar)
       
        G_n = 6.67e-8
        mstar = mstar*1.9891e33 #grams
       
 
        dim = len(sem_a)-1
        a = 0.867#1.12
        b = 3.03#3.09
        ecc = np.random.beta(a,b,size=(dim,dim,n_inject))
        print(ecc.shape)
 
        cosini = np.random.uniform(size=(dim,dim,n_inject))
        i = np.arccos(cosini)
                   
        sinii = np.sin(i)
   
        tpp = np.random.uniform(size=(dim,dim,n_inject))*(100.*365.)
 
        omm = np.random.uniform(size=(dim,dim,n_inject))*(2.*np.pi)
        ww = np.random.uniform(size=(dim,dim,n_inject))*(2.*np.pi)
        prob_s = np.zeros((dim,dim))
        prob_s_Witt = np.zeros((dim,dim))
        for l in range(len(sem_a)-1):
            print(l)
            for j in range(len(mass)-1):
                # print(j)
                tr_det = np.zeros(n_inject)
                tr_det_Witt = np.zeros(n_inject)
 
                for k in range(n_inject):
           
                    diff_a = np.log(sem_a[l+1]) - np.log(sem_a[l])
                    rand_a = np.exp((np.random.uniform(size=1))*diff_a + np.log(sem_a[l]))
                    sem_an = rand_a
                    diff_m = np.log(mass[j+1]) - np.log(mass[j])
                    rand_m = np.exp((np.random.uniform(size=1))*diff_m + np.log(mass[j]))
                    mass_n = rand_m
 
                    e = ecc[l,j,k]
                    sini = sinii[l,j,k]
                    om = omm[l,j,k]
                    w = ww[l,j,k]
                    tp = tpp[l,j,k]
 
                    sem_an = sem_an*(1.5e13) #cm
                    mass_n = mass_n*1.89e30 #grams
                    per = (4.*(np.pi**2)*(sem_an**3.)/(G_n*(mstar+(mass_n*sini))))**(0.5) #sec
                    k_amp = ((mass_n*sini)/(mstar**(2./3.)))*((per/(2.*np.pi*G_n))**(-1./3.))*(1./np.sqrt(1.-(e**2.)))/100. #m/s
                    per = per/(60.*60.*24.) #in days
                                        #print 'per'
                                        #print per
                                        #print 'kamp'
                                        #print k_amp
 
                    bigE = np.zeros(len(time))
                    nu = np.zeros(len(time))
                    for d in range(len(time)):
                        MA = 2.*(np.pi*(time[d] - tp))/per
                        bigE[d] = getE(e,MA)
                        nu[d] = 2*np.arctan((((1.+e)/(1.-e))**0.5)*np.tan(bigE[d]/2.))
                        #per = per*60.*60.*24. #sec
 
                    rv = np.zeros(len(time))
                    for h in range(len(time)):
                        rv[h] = k_amp*[(np.cos(nu[h]+w)+e*np.cos(w))]#*100. #m/s
                                                                       # print err
                    ph_err = np.random.shuffle(err)
                    ph_err = err
                    #print err
                    #print 'shuffled err'
                    #print ph_err
                    sig = np.zeros(len(time))
                    for u in range(len(time)):
                        sig[u] = np.sqrt(jitter**2 + ph_err[u])
 
                    rv_p = np.zeros(len(rv))
                    err_dist = np.zeros(len(rv))
                    for g in range(len(rv)):
                        #err_dist[g] = np.random.normal(loc=0,)*sig[g]
                        rv_p[g] = np.random.normal(loc=rv[g],scale=sig[g])#rv[g]+(err_dist[g])
 
                    #plt.figure
                    #plt.plot(time,rv_p,marker='o')
                    #plt.show()
 
                    #USE BIC
                    #single planet system
                    mod = initialize_model(per,11500,e,om,k_amp)
                    like = radvel.likelihood.RVLikelihood(mod,time,rv_p,sig)
 
                    like.params['per1'].vary = False
                    like.params['tp1'].vary = True
                    like.params['e1'].vary = False
                    like.params['w1'].vary = False
                    like.params['k1'].vary = False
                    like.params['dvdt'].vary = False
                    like.params['curv'].vary = False
                    like.params['gamma'].vary = True
                    like.params['jit'].vary = False
 
                    #print(like)
               
                    post = radvel.posterior.Posterior(like)
                    res = optimize.minimize(post.neglogprob_array,post.get_vary_params(),method='Powell')
                    ti = np.linspace(np.min(time),np.max(time),100)
       
                    res_f = post.likelihood.residuals()
        
                    chi2 = -0.5*np.sum((res_f/sig)**2.)#np.sum(((model_n - res_tr)/(err_f))**2.)
                                                   #ln_post = (-np.log(np.sqrt(2.*np.pi))*len(time_out)) -np.sum(np.log(err_out)) - (chi2[l,j,k]/(2.))
                    Likelihood_p = chi2
                    n_dat = len(time)
                    k_mod = 6.
                    BIC_p = -2.*Likelihood_p + (k_mod*np.log(n_dat))
               
 
                    #No planets, no trend
                    ones = np.ones(len(time))
                    start_par = (np.max(rv_p)+np.min(rv_p))/2.
 
                    parm = optimize.curve_fit(nopl_mod,time,rv_p,p0=start_par,sigma=sig,method='lm')
                    model_n = parm[0]*ones
                   
                    chi2_nopl = -0.5*np.sum(((rv_p - model_n)/(sig))**2.)
                    Likelihood_n = chi2_nopl
                    n_dat = len(time)
                    k_mod = 1.
                    BIC_n = -2.*Likelihood_n + (k_mod*np.log(n_dat))
 
                    if (BIC_n - BIC_p)>=10.:
                        tr_det[k] = 1.
                    else:
                        tr_det[k] = 0.
                   
                prob_s[l,j] = np.sum(tr_det)/len(tr_det)
 
        import pickle
        comp_prob = {'prob':prob_s}
        pickle.dump(comp_prob,open('complete_maps_mstar_paper/'+plname+'_complete_map.p','wb'))
       
                
def find_nearest(array,value):
    array = np.asarray(array)
    idx = (np.abs(array-value)).argmin()
    return array[idx]
 
def getE(e,M):
    Eguess = M
    deltaE = 100.
 
    while deltaE > 10.0**(-8.):
        newEguess = M + e*np.sin(Eguess)
        deltaE = np.abs(newEguess - Eguess)
        Eguess = newEguess
    return Eguess
           
def line_mod(time,p1,p2):
    rv_l = p1+time*(p2)
    return rv_l
 
def nopl_mod(time,p):
    rv_l = p+time*(0)
    return rv_l
 
def initialize_model1():
    time_base = 2420
    params = radvel.Parameters(2,basis='per tc secosw sesinw logk') # number of planets = 2
    params['per1'] = radvel.Parameter(value=20.885258)
    params['tc1'] = radvel.Parameter(value=2072.79438)
    params['secosw1'] = radvel.Parameter(value=0.01)
    params['sesinw1'] = radvel.Parameter(value=0.01)
    params['logk1'] = radvel.Parameter(value=1.1)
    params['per2'] = radvel.Parameter(value=42.363011)
    params['tc2'] = radvel.Parameter(value=2082.62516)
    params['secosw2'] = radvel.Parameter(value=0.01)
    params['sesinw2'] = radvel.Parameter(value=0.01)
    params['logk2'] = radvel.Parameter(value=1.1)
    mod = radvel.RVModel(params, time_base=time_base)
    mod.params['dvdt'] = radvel.Parameter(value=-0.02)
    mod.params['curv'] = radvel.Parameter(value=0.01)
    return mod
 
 
 
def initialize_model(per,tp,e,om,k):
    time_base = 14000
    # par,err_up,err_low = np.genfromtxt('kepler96_lm.txt',unpack=True)
 
    params = radvel.Parameters(1,basis='per tp e w k')
    #par_w = (par[3]/180.)*np.pi #convert omega to radians
    params['per1'] = radvel.Parameter(value=per)#par[0]
    params['tp1'] = radvel.Parameter(value=tp)#par[1]
    params['e1'] = radvel.Parameter(value=e)#par[2]
    params['w1'] = radvel.Parameter(value=om)#par_w
    params['k1'] = radvel.Parameter(value=k)#par[4]
 
    #print(params['per1'])
    #print(params['tp1'])
    #print(params['e1'])
    #print(params['w1'])
    #print(params['k1'])
    #print(per,tp,e,om,k)
 
    mod = radvel.RVModel(params,time_base = time_base)
    #print(mod)
    mod.params['dvdt'] = radvel.Parameter(value=0.)
    mod.params['curv'] = radvel.Parameter(value=0.)
 
    return mod
 
def plot_results(like,time,ti):
    fig = pl.figure(figsize=(12,4))
    fig = pl.gcf()
    pl.errorbar(like.x,like.model(time)+like.residuals(),yerr=like.yerr,fmt='o')
    pl.plot(ti,like.model(ti))
    fig.set_tight_layout(True)
    pl.xlabel('Time')
    pl.ylabel('RV')
    pl.draw()