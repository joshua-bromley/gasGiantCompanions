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
    
    sys = np.array(["55cnc",'upsAnd',"gamLib", "XO2S", "WASP8", "WASP53", "WASP47", "WASP41", "TYC14226141", "TOI969", "Pr0211"])#(['LSPMJ21160234'])#(['HD285968'])#(['GJ3293'])#(['GJ251'])#(['GJ229A'])#(['Wolf1061','Wolf1069','YZCet'])#(['TOI1266','TOI1452','TOI1468','TOI1695','TOI2018'])#(['proxima_cen','Ross128','Ross508','TOI244','TOI663'])#(['L9859','L36338','LHS1140','LHS1815','LTT1445A'])#(['K2-25','K2-415','Kapteyn','K2-18','Kepler138'])#(['HIP22627','HIP54373','HIP83043','HNLib','K2-3'])#(['GJ4276','GJ9689','GJ9827','HD180617','HD238090'])#(['GJ3634','GJ3779','GJ3929','GJ3988','GJ3998'])#(['GJ3082','GJ3138','GJ3323','GJ3341','GJ3470'])#(['GJ1002','GJ1132','GJ1151','GJ1214','GJ1265'])#(['GJ724','GJ740','GJ806','GJ876','GJ887'])#(['GJ674','GJ685','GJ686','GJ687','GJ720A'])#(['GJ514','GJ536','GJ581','GJ625','GJ667C'])#(['GJ411','GJ422','GJ433','GJ480','GJ486'])#(['GJ338B','GJ357','GJ367','GJ378','GJ393'])#(['GJ15A','GJ27','GJ49','GJ163','GJ273'])#(['AUMic','BD-082823','CD_Cet','G264-012','GJ12'])
   
    mass_star = np.array([0.91,1.3,1.78,0.98,1.34,0.84,1.06, 0.81, 1.15, 0.73 ])#([0.43])#([0.5])#([0.42])#([0.36])#([0.509])#([0.294,0.167,0.142])#([0.431,0.249,0.339,0.54,0.57])#([0.1221,0.168,0.1774,0.428,0.514])#([0.273,0.21,0.179,0.502,0.257])#([0.2634,0.1635,0.281,0.359,0.535])#([0.393,0.57,0.52,0.291,0.549])#([0.406,0.59,0.62,0.484,0.578])#([0.45,0.27,0.313,0.184,0.50])#([0.47,0.68,0.164,0.47,0.539])#([0.12,0.181,0.164,0.182,0.178])#([0.527,0.58,0.413,0.370,0.489])#([0.35,0.55,0.426,0.40,0.57])#([0.51,0.52,0.311,0.3,0.33])#([0.39,0.35,0.48,0.45,0.323])#([0.64,0.342,0.455,0.56,0.426])#([0.38,0.53,0.515,0.38,0.29])#([0.51,0.50,0.161,0.297,0.241])
    #mass_star = np.array([0.84])
    systemData = pd.read_csv("./RVData/rvData.txt")
    sys = systemData["hostname"].values
    mass_star = systemData["st_mass"].values

    #sys = ["BD114672"]
    #mass_star = [0.65]
    
    for x in range(305,len(sys)):
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
 
        '''
        RV, err units m/s
        Time is JD - 2450000
        '''
        '''
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
        '''
        
        #csv, Time in JD - 2400000, rv, err in m/s
        # 55 Cnc, TYC 1422-614-1
        if plname == "55cnc" or plname == "TYC14226141" or plname == "Kepler56" or plname == "HIP67851" or plname == "HD9278" or plname == "HD75898" or plname == "HD73526" or plname == "HD50499" or plname == "HD4732" or plname == "HD33844" or plname == "HD154857" or plname == "HD142" or plname == "HD137496" or plname == "HD11506" or plname == "GJ849" or plname == "GJ328" or plname == "BD202457" or plname == "gamCep" or plname == "WASP76" or plname == "TYC366712801" or plname == "TYC3318013331" or plname == "TOI2529" or plname == "TIC4672985" or plname == "HIP109600" or plname == "HIP109384" or plname == "HD89744" or plname == "HD86950" or plname == "HD5583" or plname == "HD240210" or plname == "HD238914" or plname == "HD233604" or plname == "HD222076" or plname == "HD219415" or plname == "HD19994":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time - 50000

        #csv, Time in JD, rv, err in m/s
        #ups And
        if plname == "upsAnd" or plname == "WASP47" or plname == "Pr0211" or plname == "TOI2202" or plname == 'TIC279401253' or plname == "TIC139270665" or plname == "Kepler424" or plname == "HD83443" or plname == "HD39091" or plname == "HD33142" or plname == "HD27894" or plname == "HD190360" or plname == "HD191939" or plname == "HD187123" or plname == "HD177830" or plname == "HD147873" or plname == "HD128311" or plname == "HD11964" or plname == "HD114783" or plname == "HATP13" or plname == "HATP11" or plname == "GJ317" or plname == "BD114672" or plname == "7CMa" or plname == "upsLeo" or plname == "muLeo" or plname == "gamPsc" or plname == "epsEri" or plname == "epsCrB" or plname == "betUMi" or plname == "alfTau" or plname == "alfAri" or plname == 'XO5' or plname == "WASP7" or plname == "TYC0434045831" or plname == "TOI778" or plname == "TOI6029" or plname == "TOI481" or plname == "TOI4603" or plname == "TOI2589" or plname == "TOI2497" or plname == "TOI2368" or plname== "TOI2338" or plname == "TOI2180" or plname == "TOI2107" or plname == "TOI2010" or plname == "TOI1899" or plname == "TOI1516" or plname == "TOI1181" or plname == "NGTS30":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time - 2450000
        
        if plname == "NCG2682Sand978" or plname == "NGC2682Sand1429" or plname == "Kepler43" or plname == "Kepler423" or plname == "Kepler1658" or plname == "KIC3526061" or plname == "KELT4A" or plname == "KELT24" or plname == "KELT19" or plname == "KELT17" or plname == "KELT12" or plname == "K2419A" or plname == "HSPsc" or plname == "HD99283" or plname == "HD96992" or plname == "HD9174" or plname == "HD87646" or plname == "HD81688" or plname == "HD79181" or plname == "HD68402" or plname == "HD62509" or plname == "HD52265" or plname == "HD360" or plname == "HD32963" or plname == "HD32518" or plname == "HD25723" or plname == "HD222237" or plname == "HD219139" or plname == "HD167768" or plname == "HD165155" or plname == "HD161178":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time - 2450000
        
        #csv, Time in JD, rv in km/s, err in m/s

        if  plname == "IC46519122":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time - 2450000
            rv = rv*1000

        #csv, Time in JD, rv, err in km/s
        #WASP-8, TOI-969
        if plname == "WASP8" or plname == "TOI969" or plname == "HD9446" or plname == "HD181433" or plname == "HD118203" or plname == "HD111232" or plname == "GJ676" or plname == "gam1Leo" or plname == "WASP54" or plname == "WASP42" or plname == "WASP34" or plname == "WASP3" or plname == "WASP26" or plname == "TOI5153" or plname == "NGTS20" or plname == "KIC8121913" or plname == "HD8535" or plname == "HD29021" or plname == "HD22781" or plname == "HD190984" or plname == "HD156411" or plname == "HD150706" or plname == "HD148156" or plname == "HD147513":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time - 2450000
            rv = rv*1000
            err = err*1000

        #Time in JD-2400000, rv, err in km/s
        #CoRoT24
        if plname == "CoRoT24" or plname == "Kepler539" or plname == "HIP65407" or plname == "HD219828" or plname == "HD160691" or plname == "HD141399" or plname == "HD113538" or plname == "WASP89" or plname == "WASP81" or plname == "WASP80" or plname == "WASP74" or plname == "WASP66" or plname == "WASP62" or plname == "WASP60" or plname == "WASP5" or plname == 'WASP22' or plname == "TYC4282006051" or plname == "HD42012" or plname == "HD35759" or plname == "HD331093" or plname == "HD27969" or plname == "HD233832" or plname == "HD221287" or plname == 'HD220842' or plname == "HD211403" or plname == "HD17674" or plname == "HD155193" or plname == "HD143105":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time-50000
            rv = rv*1000
            err = err*1000

        #Time in JD-2400000, rv,err in km/s, cols 0,2,3
        #XO7
        if plname == "XO7":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,2,3)).T
            time = time-50000
            rv = rv*1000
            err = err*1000

        
        #csv, Time in JD-2450000, rv,err in km/s
        #WASP-41, gam Lib
        if plname == "WASP41" or plname == "HD204313" or plname == "HD159243" or plname == "HD13908" or plname == "GJ876" or plname == "WASP73" or plname == "WASP70" or plname == "WASP68" or plname == "WASP50" or plname == "WASP4" or plname == "WASP38" or plname == "WASP32" or plname =="WASP24" or plname == "WASP106" or plname == "TOI2485" or plname == "TOI2420" or plname == "KOI3680" or plname == "CoRoT20":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            rv = rv*1000
            err = err*1000


        #csv, Time in JD-2450000, rv,err in m/s
        if plname == "gamLib" or plname == "WASP48" or plname == "TrES4" or plname == "TrES2" or plname == "Qatar2" or plname == "HIP107773" or plname == "HD76920" or plname == "HD66141" or plname == "HD60292" or plname == "HD59686A" or plname == "HD47526" or plname == "HD36384" or plname == "HD2952" or plname == "HD220773" or plname == "HD210702" or plname == "HD208897" or plname == "HD208527" or plname == "HD197037" or plname == "HD19615" or plname == "HD174205" or plname == "HD17416" or plname == "HD167402" or plname == "HD1666" or plname == "HD154391" or plname == "HD150010" or plname == "HD145457":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
        
        if plname == "Kepler93" or plname == "Kepler94" or plname == "Kepler68" or plname == "Kepler48" or plname == "Kepler454" or plname == "Kepler407" or plname == "KELT6" or plname == "HIP8541" or plname == "HD80653" or plname == "HD67087" or plname == "HD47366" or plname == "HD155358" or plname == "HD134987" or plname == "HATS59" or plname == "HATP2" or plname == "75Cet" or plname == "xiAql" or plname == "tauBoo" or plname == "omeSer" or plname == "HIP75092" or plname == "HIP74890" or plname == "HIP114933" or plname == "HD158996" or plname == "HD14067":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T

         #time in BJD-2457000, rv, err in m/s
        if plname == "TOI1130":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time + 7000   
        #time in BJD-2455000, rv,err in m/s
        if plname == "Kepler432" or plname == "HIP105854":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time + 5000
        
        #time in JD-2440000, rv, err in m/s
        if plname == "HIP14810" or plname == "HD99706" or plname == "HD75784" or plname == "HD74156" or plname == "HD5319" or plname == "HD4203" or plname == "HD37605" or plname == "HD37124" or plname == "HD30177" or plname == "HD207832" or plname == "HD200964" or plname == "HD1605" or plname == "HD159868" or plname == "HD148164" or plname == "HD125612" or plname == "HD116029" or plname == "HD108874" or plname == "HD102329" or plname == "HATP17" or plname == "47UMa" or plname == "24Sex" or plname == "HR5183" or plname == "HD98736" or plname == "HD98219" or plname == "HD96167" or plname == "HD96063" or plname == "HD94834" or plname == "HD86081" or plname == "HD73534" or plname == "HD72490" or plname == "HD55696" or plname == "HD4917":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time - 10000
        if plname == "HD43691" or plname == "HD4313" or plname == "HD30856" or plname == "HD28678" or plname == "HD231701" or plname == "HD224693" or plname == "HD212771" or plname == "HD211810" or plname == "HD206610" or plname == "HD196885A" or plname == "HD18742" or plname == "HD185269" or plname == "HD181342" or plname == "HD180902" or plname == "HD18015" or plname == "HD180053" or plname == "HD175541" or plname == "HD17156" or plname == "HD152581" or plname == "HD1502" or plname == "HD149143" or plname == "HD14787":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time - 10000

        #time in JD, rv,err in m/s, CLS tables
        if plname == "HD72659" or plname == "HD68988" or plname == "HD38529" or plname == "HD28185" or plname == "HD24040" or plname == "HD183263" or plname == "HD169830" or plname == "HD168443" or plname == "HD156279" or plname == "HD12661" or plname == "14Her" or plname == "HD87883" or plname == "HD8574" or plname == "HD80606" or plname == "HD50554" or plname == "HD213472" or plname == "HD209458" or plname == "HD189733" or plname == "HD136925":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (5,6,7)).T
            time = time - 2450000

        #time in JD-2454000, rv, err in km/s
        if plname == "HD113337":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time + 4000
            rv = rv*1000
            err = err*1000

        #time in JD-2453000, rv, err in km/s
        if plname == "HD221585":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time + 3000
            rv = rv*1000
            err = err*1000

        #time in JD-2453000, rv, err in m/s
        if plname == "HD205739":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time + 3000


        #time in JD-2456000, rv, err in km/s
        if plname == "Kepler447":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time + 6000
            rv = rv*1000
            err = err*1000

        #time in JD-2454000, rv, err in m/s
        if plname == "HATP44":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time + 4000
        #time in JD-2454900, rv,err in m/s
        if plname == "KOI142":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
            time = time + 4900

        #tsv, time in JD-2450000, rv err in km/s
        #XO-2 S
        if plname == "XO2S" or plname == "CoRoT20":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter="\t", usecols = (0,1,2)).T
            rv = rv*1000
            err = err*1000

        #tsv, time in JD-2400000, rv, err in km/s
        #WASP 53
        if plname == "WASP53":
            time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter="\t", usecols = (0,1,2)).T
            time = time - 50000
            rv = rv*1000
            err = err*1000
 
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
        pickle.dump(comp_prob,open(f'./completenessMaps/{plname}.p','wb'))
       
                
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

complete_maps()