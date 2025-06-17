import pandas as pd
import numpy as np

planets = pd.read_csv("./data/gasGiantDataComplete.csv")

baseline = np.zeros(len(planets))
nObs = np.zeros(len(planets))

planets["baseline"] = baseline
planets["num_obs"] = nObs

for i in range(len(planets)):
    planet = planets.iloc[i]
    plname = planet["hostname"]

    if plname == "55Cnc" or plname == "TYC14226141" or plname == "Kepler56" or plname == "HIP67851" or plname == "HD9278" or plname == "HD75898" or plname == "HD73526" or plname == "HD50499" or plname == "HD33844" or plname == "HD154857" or plname == "HD142" or plname == "HD137496" or plname == "HD11506" or plname == "GJ849" or plname == "GJ328" or plname == "BD202457" or plname == "gamCep" or plname == "WASP76" or plname == "TYC366712801" or plname == "TYC3318013331" or plname == "TOI2529" or plname == "TIC4672985" or plname == "HIP109600" or plname == "HIP109384" or plname == "HD89744" or plname == "HD86950" or plname == "HD5583" or plname == "HD240210" or plname == "HD238914" or plname == "HD233604" or plname == "HD222076" or plname == "HD219415" or plname == "HD19994" or plname == "BD490828" or plname == "BD480740" or plname == 'BD480738' or plname == "BD200274" or plname == "BD152940" or plname == "BD152375" or plname == "BD032562":
        time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
        time = time - 50000

        #csv, Time in JD, rv, err in m/s
        #ups And
    if plname == "upsAnd" or plname == "WASP47" or plname == "Pr0211" or plname == "TOI2202" or plname == 'TIC279401253' or plname == "TIC139270665" or plname == "Kepler424" or plname == "HD83443" or plname == "HD39091" or plname == "HD33142" or plname == "HD27894" or plname == "HD190360" or plname == "HD191939" or plname == "HD187123" or plname == "HD177830" or plname == "HD147873" or plname == "HD128311" or plname == "HD11964" or plname == "HD114783" or plname == "HATP13" or plname == "HATP11" or plname == "GJ317" or plname == "BD114672" or plname == "7CMa" or plname == "upsLeo" or plname == "muLeo" or plname == "gamPsc" or plname == "epsEri" or plname == "epsCrB" or plname == "betUMi" or plname == "alfTau" or plname == "alfAri" or  plname == "WASP7" or plname == "TYC0434045831" or plname == "TOI778" or plname == "TOI6029" or plname == "TOI481" or plname == "TOI4603" or plname == "TOI2589" or plname == "TOI2497" or plname == "TOI2368" or plname== "TOI2338" or plname == "TOI2180" or plname == "TOI2107" or plname == "TOI2010" or plname == "TOI1899" or plname == "TOI1516" or plname == "TOI1181" or plname == "NGTS30" or plname == "xiAql" or plname == "omeSer" or plname == "gamLib" or plname == "HD210702" or plname == "HD208897" or plname == "HD167042" or plname == "HD14067" or plname == "HD104985"  or plname == "6Lyn" or plname == "HATP17":
        time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
        time = time - 2450000
        
    if plname == "NCG2682Sand978" or plname == "NGC2682Sand1429" or plname == "Kepler43" or plname == "Kepler423" or plname == "Kepler1658" or plname == "KIC3526061" or plname == "KELT4A" or plname == "KELT24" or plname == "KELT19" or plname == "KELT17" or plname == "KELT12" or plname == "K2419A" or plname == "HSPsc" or plname == "HD99283" or plname == "HD96992" or plname == "HD9174" or plname == "HD87646" or plname == "HD81688" or plname == "HD79181" or plname == "HD68402" or plname == "HD62509" or plname == "HD52265" or plname == "HD360" or plname == "HD32963" or plname == "HD32518" or plname == "HD25723" or plname == "HD222237" or plname == "HD219139" or plname == "HD167768" or plname == "HD165155" or plname == "HD161178" or plname == "HD132563" or plname == "HD128356" or plname == "HD114082" or plname == "HATS70" or plname == "HATS56" or plname == "HATS41" or plname == "HATP34" or plname == "HATP33" or plname == "HATP32" or plname == "HATP31" or plname == "HATP30" or plname == "HATP29" or plname == "HATP24" or plname == "HATP16" or plname == "HATP15" or plname == "4UMa" or plname == "24Boo" or plname == "18Del" or plname == "17Sco" or plname == "14And" or plname == "11Com" or plname == "HD47366" or plname == "HD4732" or plname == "HD2952" or plname == "HD173416" or plname == "HD145457" or plname == "HD120084" or plname == "81Cet":
        time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
        time = time - 2450000
        
        #csv, Time in JD, rv in km/s, err in m/s

    if  plname == "IC46519122":
        time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
        time = time - 2450000
        rv = rv*1000

        #csv, Time in JD, rv, err in km/s
        #WASP-8, TOI-969
    if plname == "WASP8" or plname == "TOI969" or plname == "HD9446" or plname == "HD181433" or plname == "HD118203" or plname == "HD111232" or plname == "GJ676" or plname == "gam1Leo" or plname == "WASP54" or plname == "WASP42" or plname == "WASP34" or plname == "WASP3" or plname == "WASP26" or plname == "TOI5153" or plname == "NGTS20" or plname == "KIC8121913" or plname == "HD8535" or plname == "HD29021" or plname == "HD22781" or plname == "HD190984" or plname == "HD156411" or plname == "HD150706" or plname == "HD148156" or plname == "HD147513" or plname == "HATP70" or plname == "HATP69" or plname == "GJ3021" or plname == "CoRoT29" or plname == "CoRoT12" or plname == "42Dra":
        time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
        time = time - 2450000
        rv = rv*1000
        err = err*1000

        #Time in JD-2400000, rv, err in km/s
        #CoRoT24
    if plname == "CoRoT24" or plname == "Kepler539" or plname == "HIP65407" or plname == "HD219828" or plname == "HD160691" or plname == "HD141399" or plname == "HD113538" or plname == "WASP89" or plname == "WASP81" or plname == "WASP80" or plname == "WASP74" or plname == "WASP66" or plname == "WASP62" or plname == "WASP60" or plname == "WASP5" or plname == 'WASP22' or plname == "TYC4282006051" or plname == "HD42012" or plname == "HD35759" or plname == "HD331093" or plname == "HD27969" or plname == "HD233832" or plname == "HD221287" or plname == 'HD220842' or plname == "HD211403" or plname == "HD17674" or plname == "HD155193" or plname == "HD143105" or plname == "HD12484" or plname == "HD124330" or plname == "HD115954" or plname == "HD109286" or plname == 'HD108341' or plname == "HD103720" or plname == "CoRoT26" or plname == "BD550362":
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
    if plname == "WASP41" or plname == "HD204313" or plname == "HD159243" or plname == "HD13908" or plname == "GJ876" or plname == "WASP73" or plname == "WASP70" or plname == "WASP68" or plname == "WASP50" or plname == "WASP4" or plname == "WASP38" or plname == "WASP32" or plname =="WASP24" or plname == "WASP106" or plname == "TOI2485" or plname == "TOI2420" or plname == "KOI3680" or plname == "CoRoT20" or plname == "HATS54" or plname == "HATP3" or plname == "HATP21" or plname == "HATP20" or plname == "HATP14" or plname == "CoRoT9":
        time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
        rv = rv*1000
        err = err*1000


        #csv, Time in JD-2450000, rv,err in m/s
    if plname == "WASP48" or plname == "TrES4" or plname == "TrES2" or plname == "Qatar2" or plname == "HIP107773" or plname == "HD76920" or plname == "HD66141" or plname == "HD60292" or plname == "HD59686A" or plname == "HD47526" or plname == "HD36384" or plname == "HD220773" or plname == "HD208527" or plname == "HD197037" or plname == "HD19615" or plname == "HD174205" or plname == "HD17416" or plname == "HD167402" or plname == "HD1666" or plname == "HD154391" or plname == "HD150010" or plname == "HD11977" or plname == "HD113996" or plname == "HD112640" or plname == "HD112570" or plname == "HD112300" or plname == "HATS51" or plname == "HATS50" or plname == "HATS33":
        time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
        
    if plname == "Kepler93" or plname == "Kepler94" or plname == "Kepler68" or plname == "Kepler48" or plname == "Kepler454" or plname == "Kepler407" or plname == "KELT6" or plname == "HIP8541" or plname == "HD80653" or plname == "HD67087" or plname == "HD155358" or plname == "HD134987" or plname == "HATS59" or plname == "HATP2" or plname == "75Cet" or plname == "tauBoo" or plname == "HIP75092" or plname == "HIP74890" or plname == "HIP114933" or plname == "HD158996" or plname == "HD110014" or plname == "HATS68" or plname == "HATS64" or plname == "HATS62" or plname == "HATS27" or plname == "HATP66" or plname == "HATP63" or plname == "HATP61" or plname == "HATP60" or plname == "HATP22" or plname == "HATP1" or plname == "91Aqr":
        time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T

         #time in BJD-2457000, rv, err in m/s
    if plname == "TOI1130":
        time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
        time = time + 7000   
        #time in BJD-2455000, rv,err in m/s
    if plname == "Kepler432" or plname == "HIP105854":
        time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
        time = time + 5000
        #time in BJD-2455000, rv,err in km/s
    if plname == "CoRoT16":
        time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
        time = time + 5000
        rv = rv*1000
        err = err*1000
        
        #time in JD-2440000, rv, err in m/s
    if plname == "HIP14810" or plname == "HD99706" or plname == "HD75784" or plname == "HD74156" or plname == "HD5319" or plname == "HD4203" or plname == "HD37605" or plname == "HD37124" or plname == "HD30177" or plname == "HD207832" or plname == "HD200964" or plname == "HD1605" or plname == "HD159868" or plname == "HD148164" or plname == "HD125612" or plname == "HD116029" or plname == "HD108874" or plname == "HD102329"  or plname == "47UMa" or plname == "24Sex" or plname == "HR5183" or plname == "HD98736" or plname == "HD98219" or plname == "HD96167" or plname == "HD96063" or plname == "HD94834" or plname == "HD86081" or plname == "HD73534" or plname == "HD72490" or plname == "HD55696" or plname == "HD4917":
        time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
        time = time - 10000
    if plname == "HD43691" or plname == "HD4313" or plname == "HD30856" or plname == "HD28678" or plname == "HD231701" or plname == "HD224693" or plname == "HD212771" or plname == "HD211810" or plname == "HD206610" or plname == "HD196885A" or plname == "HD18742" or plname == "HD185269" or plname == "HD181342" or plname == "HD180902" or plname == "HD18015" or plname == "HD180053" or plname == "HD175541" or plname == "HD17156" or plname == "HD152581" or plname == "HD1502" or plname == "HD149143" or plname == "HD14787" or plname == "HD131496" or plname == "HD130322" or plname == "HD108863" or plname == "HD106270" or plname == "HD10442" or plname == "HD102956" or plname == "GJ179" or plname == "70Vir" or plname == 'XO5':
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
    if plname == "HATP44" or plname == "HATP50" or plname == "HATP39":
        time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
        time = time + 4000
        #time in JD-2454900, rv,err in m/s
    if plname == "KOI142":
        time, rv, err = np.loadtxt(f"./RVData/rv{plname}.txt", delimiter=",", usecols = (0,1,2)).T
        time = time + 4900

        #tsv, time in JD-2450000, rv err in km/s
        #XO-2 S
    if plname == "XO2S":
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
    try:
        baseline[i] = np.max(time) - np.min(time)
        nObs[i] = len(time)
    except NameError:
        print(plname)

planets["baseline"] = baseline
planets["num_obs"] = nObs


planets.to_csv("./data/gasGiantDataComplete.csv")    
