import numpy as np

class SlabOut:
    def __init__(self, dir_out, spec_list=None, 
            NH_total=None, Zd=1., xC=1.6e-4, xO=3.2e-4, xSi=1.7e-6, xS=3.5e-6):
        self.dir_out = dir_out
        if spec_list == None:
            self.spec_list = ["He+", "OHx", "CHx", "CO",  "C+", "HCO+", "H2", "H+", "H3+", "H2+","S+","Si+", "O+", "E"]
        else:
            self.spec_list = spec_list
        self.xC = xC
        self.xO = xO
        print("Set xC={:.2e}, xO={:.2e}.".format(self.xC, self.xO))
        self.nH = np.loadtxt(dir_out+"nH_arr.dat")
        self.nslab = len(self.nH)
        slab_arr = []
        for islab in np.arange(self.nslab):
            fn_slab = dir_out + "slab" + ("%06d" % islab) + ".dat"
            slab_arr.append(np.loadtxt(fn_slab))
        slab_arr = np.array(slab_arr)
        self.abd = {}
        for i in np.arange(len(self.spec_list)):
            self.abd[self.spec_list[i]] = np.transpose(slab_arr[:, :, i])
        self.abd["C"] = self.xC - self.abd["CHx"] - self.abd["CO"] - self.abd["C+"] - self.abd["HCO+"]
        self.abd["2H2"] = self.abd["H2"]*2
        self.abd["H"] = 1. - (self.abd["H2"]*2 + 3*self.abd["H3+"] + 2*self.abd["H2+"] 
                              + self.abd["H+"] + self.abd["HCO+"] + self.abd["CHx"] + self.abd["OHx"])
        self.abd["CO2Ctot"] = self.abd["CO"]/self.xC
        self.abd["C+2Ctot"] = self.abd["C+"]/self.xC
        self.abd["C2Ctot"] = self.abd["C"]/self.xC
        self.abd["e"] = self.abd["He+"] + self.abd["C+"] + self.abd["HCO+"] + self.abd["H+"
                        ] + self.abd["H3+"] + self.abd["H2+"]
        self.abd["O"] = self.xO - self.abd["OHx"] - self.abd["CO"] - self.abd["HCO+"]
        self.abd["He"] = 0.1 - self.abd["He+"]
        self.abd["M+"] = np.zeros(self.abd["C"].shape)
        if "S+" in self.spec_list:
            self.abd["M+"] = self.abd["M+"] + self.abd["S+"]
            self.abd["e"] = self.abd["e"] + self.abd["S+"]
        if "Si+" in self.spec_list:
            self.abd["M+"] = self.abd["M+"] + self.abd["Si+"]
            self.abd["e"] = self.abd["e"] + self.abd["Si+"]
        if xS != None:
            self.abd["S"] = xS - self.abd["S+"]
        if xSi != None:
            self.abd["Si"] = xSi - self.abd["Si+"]
        if "E" in self.spec_list:
            print("Calculating E assuming CvCold and xHe=0.1 ...")
            kb = 1.381e-16
            self.abd["T"] = self.abd["E"] / (1.5 * kb * (self.abd["H2"] + 1. - 2.*self.abd["H2"] + self.abd["e"] + 0.1))
        self.ngrid = slab_arr.shape[1]
        if NH_total != None:
            self.NH_total = NH_total
            self.NH = np.linspace(0, self.NH_total, self.ngrid)
        else:
            self.NH = np.loadtxt(dir_out+"colH_arr.dat")
        self.Zd = Zd
        self.Av = self.NH * Zd/1.87e21
        self.nHM, self.NHM = np.meshgrid(self.nH, self.NH)
        self.NH2 = np.zeros( (self.ngrid, self.nslab) )
        self.NCO = np.zeros( (self.ngrid, self.nslab) )
    def GetAbd(self, spec, nH=None, NH=None):
        if nH==None and NH==None:
            return self.abd[spec]
        elif nH!=None and NH==None:
            indx_nH = np.argmin(abs(self.nH - nH))
            return self.abd[spec][:, indx_nH]
        elif nH==None and NH!=None:
            indx_NH = np.argmin(abs(self.NH - NH))
            return self.abd[spec][indx_NH,:]
        else:
            indx_nH = np.argmin(abs(self.nH - nH))
            indx_NH = np.argmin(abs(self.NH - NH))
            return self.abd[spec][indx_NH, indx_nH]
    def ReadThermo(self, thermo_list=None):
        if thermo_list == None:
            self.thermo_list = ["LCR" , "LPE" , "LH2gr" , "LH2pump" , "LH2diss",
                                "GCII" , "GCI" , "GOI" , "GLya" , "GCOR",
                                "GH2" , "GDust" , "GRec" , "GH2diss" , "GHIion"]
        else:
            self.thermo_list = thermo_list
        thermo_arr = []
        for islab in np.arange(self.nslab):
            fn_thermo = self.dir_out + "thermo" + ("%06d" % islab) + ".dat"
            thermo_arr.append(np.loadtxt(fn_thermo))
        thermo_arr = np.array(thermo_arr)
        for i in np.arange(len(self.thermo_list)):
            self.abd[self.thermo_list[i]] = np.transpose(thermo_arr[:, :, i])
        self.abd["Ltotal"] = np.zeros(self.abd["T"].shape)
        self.abd["Gtotal"] = np.zeros(self.abd["T"].shape)
        for x in self.thermo_list:
            if x[0] == "L":
                self.abd["Ltotal"] += self.abd[x]
            elif x[0] == "G":
                self.abd["Gtotal"] += self.abd[x]
            else:
                print("ERROR: can't identify process {}".format(x))
    def ReadRates(self, rates_list=None):
        if rates_list == None:
            self.ratefile = self.dir_out + "chemnet.dat"
            f = open(self.ratefile)
            flines = f.readlines()
            self.rates_list = [fl.split(",")[0] for fl in flines]
        else:
            self.rates_list = rates_list
        rates_arr = []
        for islab in np.arange(self.nslab):
            fn_rates = self.dir_out + "rates" + ("%06d" % islab) + ".dat"
            rates_arr.append(np.loadtxt(fn_rates))
        rates_arr = np.array(rates_arr)
        self.rates = {}
        self.nrates = len(self.rates_list)
        for i in np.arange(self.nrates):
            self.rates[i] = np.transpose(rates_arr[:, :, i])
            self.abd[self.rates_list[i]] = np.transpose(rates_arr[:, :, i])         


