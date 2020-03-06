import numpy as np
import matplotlib.pyplot as plt

#reactions for one species formation/destruction
def GetReaction(reaction_list, spec):
    creation = []
    destruction = []
    for r in reaction_list:
        rsplit = r.split(" ")
        if spec in rsplit:
            s = r.split("->")
            sp0 = s[0].split(" ")
            sp1 = s[1].split(" ")
            if spec in sp0:
                destruction.append(r)
            else:
                creation.append(r)
    return creation, destruction

def GetNi(NH_arr, xi_arr):
    if NH_arr.shape != xi_arr.shape:
        raise Exception("GetNi: array shape doesn't match.")
    Ni = np.zeros(NH_arr.shape)
    Ni[0]  = xi_arr[0] * NH_arr[0]
    for i in np.arange(1, len(Ni)):
        Ni[i] = Ni[i-1] + xi_arr[i] * (NH_arr[i] - NH_arr[i-1])
    return Ni

def checkabd(slab):
    Ctot = (slab.abd["C"] + slab.abd["C+"] + slab.abd["CO"])[0,0]
    Otot = (slab.abd["O"] + slab.abd["CO"])[0,0]
    Stot = (slab.abd["S"] + slab.abd["S+"])[0,0]
    Sitot = (slab.abd["Si"] + slab.abd["Si+"])[0,0]
    PAHtot = slab.abd["PAHtotal"][0, 0]
    print("Ctot={:.2e}, Otot={:.2e}, Stot={:.2e}, Sitot={:.2e}, PAHtot={:.2e}".format(
        Ctot, Otot, Stot, Sitot, PAHtot))
    return

def plot_Ntrans_BS16(ax):
    n = np.logspace(0, 4, 100)
    RFUV = 5.7e-11
    sigma_g = 2.24e-21
    sigma_g_tilt = sigma_g / 1.9e-21
    pre_factor = 0.59 * (RFUV/5.8e-11) * (9.9/(1. + 8.9*sigma_g_tilt))
    alphaG = pre_factor * (100./n)
    N = 0.7/sigma_g * np.log( (alphaG/2.)**1.43 + 1 )
    line=ax.plot(N, n, "k-.")
    line[0].set_dashes([10, 4, 2.2, 4])
    return

def Get_xOH_BS15(xi_H, n, x_O, x_e, Z):
    xi_2 = xi_H * 2.3
    return (xi_2/n/1e-16) * (1 + 0.27*x_O/x_e) / Z

def plot_contour(plotslabs, linestyles=["-", "--", ":"], figname_end="", savefig=True, plot_BS16=False):
    plot_spec_list=["CO2Ctot", "C+2Ctot", "2H2"]
    plot_label=["$x(\mathrm{CO})/x(\mathrm{C_{tot}})$", "$x(\mathrm{C^+})/x(\mathrm{C_{tot}})$", 
                "$2x(\mathrm{H_2})$"]
    colors = ["m", "y", "k"]
    fig_dir = "/Users/munangong/chemistry_Athena/draft_chemistry/"
    figname="contour" + figname_end
    level = 0.5
    fig=plt.figure(figsize=[10,8])
    ax = fig.add_subplot(111)
    for islab in np.arange(len(plotslabs)):
        for i in np.arange(len(plot_spec_list)):
            CS1 = ax.contour(plotslabs[islab].NHM, plotslabs[islab].nHM, plotslabs[islab].abd[plot_spec_list[i]], 
                             [level], colors=colors[i], linestyles=linestyles[islab])
            if linestyles[islab] == "--":
                CS1.collections[0].set_dashes([(0, (20.0, 8.0))])
            if linestyles[islab] == ":":
                CS1.collections[0].set_dashes([(0, (4.2, 6.0))])

            if islab == 0:
                CS1.collections[0].set_label(plot_label[i])
    high_line = np.ones(plotslabs[islab].NH.shape) * 1e4
    ax.fill_between(plotslabs[islab].NH, 31./plotslabs[islab].NH * 100 * 1e21, high_line, facecolor='gray', alpha=0.5)
    ax.fill_between(plotslabs[islab].NH, 3.1/plotslabs[islab].NH * 100 * 1e21, high_line, facecolor='gray', alpha=0.2)
    if plot_BS16:
        plot_Ntrans_BS16(ax)
    #ax.plot(plotslabs[islab].NH, 31./plotslabs[islab].NH * 100 * 1e21, "b-")
    #ax.plot(plotslabs[islab].NH, 3.1/plotslabs[islab].NH * 100 * 1e21, "b--")
    #ax.text(3.5e20, 145, r"$t_\mathrm{H_2}/t_\mathrm{dyn}<10$", fontsize=25)
    #ax.text(1.4e21, 449, r"$t_\mathrm{H_2}/t_\mathrm{dyn}<1$", fontsize=25)

    ax.set_ylabel("$n/\mathrm{cm^{-3}}$", fontsize=25)
    ax.set_xlabel("$N/\mathrm{cm^{-2}}$", fontsize=25)
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylim([50, 1000])
    ax.set_xlim([1e19, 1e22])
    ax.legend(loc=2, fontsize=23)

    ax3 = ax.twiny()
    ax3.set_xlim(np.array(ax.get_xlim())/1.87e21)
    ax3.set_xscale(ax.get_xscale())
    ax3.set_xlabel("$A_V$", fontsize=25)
    fig.tight_layout()
    if savefig:
        fig.savefig(fig_dir + figname + ".pdf")
    return


def plot_xi_simple(plot_slabs, linestyles=["-", "--", ":"], figname_end="", savefig=True, plot_nH_arr=[100, 1000], Z=1.,
                   plot_xOH_BS15=False, xi_H_arr=None):
    #plot in results, aboundances, cr rate 2e-16
    fig_dir = "/Users/munangong/chemistry_Athena/draft_chemistry/"

    if "TTB" in figname_end:
        plot_spec_list=["CO", "C", "C+", "H3+", "OHx", "CHx", "He+", "T"]
    else:
        plot_spec_list=["CO", "C", "C+", "H3+", "OHx", "CHx", "He+"]
    plot_spec_label=["$\mathrm{CO}$", "$\mathrm{C}$", "$\mathrm{C^+}$", "$\mathrm{H_3^+}$", 
                     "$\mathrm{OH_x}$", "$\mathrm{CH_x}$", "$\mathrm{He^+}$", r"$10^{-6}\mathrm{T}$"]
    colors = ["m", "r", "y", "g", "b", "cyan", "k", "gray"]
    indxs = np.arange(len(plot_spec_list))

    for inH in np.arange(len(plot_nH_arr)):
        plot_nH = int(plot_nH_arr[inH])
        figname = "species" +figname_end + "_nH"+str(plot_nH)
        fig=plt.figure(figsize=[10,11])
        ax = fig.add_subplot(111)
        for islab in np.arange(len(plot_slabs)):
            for s,i in zip(plot_spec_list, indxs):
                if islab == 0:
                    slabel = plot_spec_label[i]
                else:
                    slabel = None
                if s == "T":
                    factor=1e-6
                else:
                    factor=1
                line = ax.plot(plot_slabs[islab].NH, plot_slabs[islab].GetAbd(s, nH=plot_nH)*factor, label=slabel, 
                               color=colors[i], linestyle=linestyles[islab])
                if linestyles[islab] == "-.":
                    line[0].set_dashes([20, 6, 5, 6])
                if linestyles[islab] == "--":
                    line[0].set_dashes([20, 8, 20, 8])
                if linestyles[islab] == ":":
                    line[0].set_dashes([4.2, 6, 4.2, 6])

        ax.set_xlim([5e19/Z, 1e22/Z]);
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend(fontsize=20, loc=3)
        ax.set_xlabel("$N/\mathrm{cm^{-2}}$", fontsize=25)
        ax.set_ylabel("$x_s = n_s/n$", fontsize=25)
        ax3 = ax.twiny()
        ax3.set_xlim(np.array(ax.get_xlim())/1.87e21*Z)
        ax3.set_xscale(ax.get_xscale())
        ax3.set_xlabel("$A_V$", fontsize=25)
        ax3.set_ylim([1e-9, 2e-4])
        ax.xaxis.set_ticks_position('bottom')
        #ax.set_title("nH={:d}\n\n\n\n".format(plot_nH), fontsize=18)
        if plot_xOH_BS15:
            x_O = plot_slabs[islab].GetAbd("O", nH=plot_nH, NH=1e22/Z)
            x_e = plot_slabs[islab].GetAbd("e", nH=plot_nH, NH=1e22/Z)
            xOH_BS15 = Get_xOH_BS15(xi_H_arr[inH], float(plot_nH), x_O, x_e, Z)
            ax.plot(1e22/Z, xOH_BS15, "b*", ms=10, linestyle="None")
            print("nH={}, x_O={}, x_e={}, xOH_BS15={}".format(plot_nH, x_O,
                x_e, xOH_BS15))
        fig.tight_layout()
        if savefig:
            fig.savefig(fig_dir + figname + ".pdf")
    return
            
def plot_Ni_simple(plot_slabs, linestyles=["-", "--", ":"], figname_end="", savefig=True, plot_nH_arr=[100, 1000]):           
    fig_dir = "/Users/munangong/chemistry_Athena/draft_chemistry/"

    plot_spec_list=["CO", "C", "C+", "H3+", "OHx", "CHx", "He+"]
    plot_spec_label=["$\mathrm{CO}$", "$\mathrm{C}$", "$\mathrm{C^+}$", "$\mathrm{H_3^+}$", 
                     "$\mathrm{OH_x}$", "$\mathrm{CH_x}$", "$\mathrm{He^+}$"]
    colors = ["m", "r", "y", "g", "b", "cyan", "k"]

    for plot_nH in plot_nH_arr:
        plot_nH = int(plot_nH)
        figname = "Ni" +figname_end + "_nH"+str(plot_nH)
        fig=plt.figure(figsize=[10,11])
        ax = fig.add_subplot(111)
        for islab in np.arange(len(plot_slabs)):
            for s,i in zip(plot_spec_list, np.arange(len(plot_spec_list))):
                if islab == 0:
                    slabel = plot_spec_label[i]
                else:
                    slabel = None
                line = ax.plot(plot_slabs[islab].NH, 
                               GetNi(plot_slabs[islab].NH, plot_slabs[islab].GetAbd(s, nH=plot_nH)), label=slabel, 
                               color=colors[i], linestyle=linestyles[islab])
                if linestyles[islab] == "-.":
                    line[0].set_dashes([10, 4, 2.2, 4])
                if linestyles[islab] == "--":
                    line[0].set_dashes([10, 5, 10, 5])
                if linestyles[islab] == ":":
                    line[0].set_dashes([2.2, 4, 2.2, 4])

        ax.set_xlim([5e19, 1e22]);
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend(fontsize=20, loc=2)
        ax.set_xlabel("$N/\mathrm{cm^{-2}}$", fontsize=25)
        ax.set_ylabel("$N_i/\mathrm{cm^{-2}}$", fontsize=25)
        ax3 = ax.twiny()
        #ax3.plot(plot_slab1.NH/1.87e21, plot_slab1.GetAbd("e", nH=plot_nH), color="gray", linestyle="none")
        ax3.set_xlim(np.array(ax.get_xlim())/1.87e21)
        ax3.set_xscale(ax.get_xscale())
        ax3.set_xlabel("$A_V$", fontsize=25)
        ax3.set_ylim([1e12, 2e18])
        ax.xaxis.set_ticks_position('bottom')
        #ax.set_title("nH={:d}\n\n\n\n".format(plot_nH), fontsize=18)
        fig.tight_layout()
        if savefig:
            fig.savefig(fig_dir + figname + ".pdf")
    return
            
#all species in appendix
def Plot_xi(plot_slabs, plot_nH_arr, Zd, figname_end="", savefig=False, linestyles=["-", "--", ":"]):
    fig_dir = "/Users/munangong/chemistry_Athena/draft_chemistry/"

    if "OH" in figname_end:
        plot_spec_list=["CO", "C", "C+", "HCO+", "OHx", "CHx", "O", "O+", "OH+", "H2O+"]
    else:
        plot_spec_list=["CO", "C", "C+", "HCO+", "OHx", "CHx", "O", "O+"]
    plot_spec_label=["$\mathrm{CO}$", "$\mathrm{C}$", "$\mathrm{C^+}$", "$\mathrm{HCO^+}$", 
                     "$\mathrm{OH_x}$", "$\mathrm{CH_x}$", "$\mathrm{O}$", 
                     "$\mathrm{O^+}$", "$\mathrm{OH^+}$", "$\mathrm{H_2O^+}$"]
    colors = ["m", "r", "y", "g", "b", "cyan", "k", "gray", "BurlyWood", "DarkCyan"]

    plot_spec_list2=["2H2", "H+", "H2+", "H3+", "He+", "Si+", "e"]
    if "TTB" in figname_end:
        plot_spec_list2=["2H2", "H+", "H2+", "H3+", "He+", "Si+", "e", "T"]
    plot_spec_label2=[r"$2\times 10^{-3}\mathrm{H_2}$", "$\mathrm{H^+}$", "$\mathrm{H_2^+}$",
                      "$\mathrm{H_3^+}$", "$\mathrm{He^+}$",
                      "$\mathrm{Si^+}$", "$\mathrm{e}$", r"$10^{-6}\mathrm{T}$"]
    colors2 = ["k", "y", "g", "b", "r", "cyan", "m", "gray"]

    for plot_nH in plot_nH_arr:
        plot_nH = int(plot_nH)
        figname = "species_all_nH"+str(plot_nH) + "Z_" + str(int(Zd)) + figname_end
        fig=plt.figure(figsize=[20,8])
        ax = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        for islab in np.arange(len(plot_slabs)):
            for s,i in zip(plot_spec_list, np.arange(len(plot_spec_list))):
                if islab == 0:
                    slabel = plot_spec_label[i]
                else:
                    slabel = None
                factor = 1.
                try:
                    line = ax.plot(plot_slabs[islab].NH, plot_slabs[islab].GetAbd(s, nH=plot_nH)*factor, label=slabel, 
                                   color=colors[i], linestyle=linestyles[islab])
                    if linestyles[islab] == "-.":
                        line[0].set_dashes([10, 4, 2.2, 4])
                    if linestyles[islab] == "--":
                        line[0].set_dashes([10, 5, 10, 5])
                    if linestyles[islab] == ":":
                        line[0].set_dashes([2.5, 4, 2.5, 4])
                except KeyError:
                    pass
                if s == "CHx" or s == "OHx":
                    try:
                        line = ax.plot(plot_slabs[islab].NH, 
                                       plot_slabs[islab].GetAbd(s[:-1], nH=plot_nH), 
                                       color=colors[i], linestyle="-.")
                        line[0].set_dashes([10, 4, 2.2, 4])
                    except KeyError:
                        pass
            for s,i in zip(plot_spec_list2, np.arange(len(plot_spec_list2))):
                if islab == 0:
                    slabel = plot_spec_label2[i]
                else:
                    slabel = None
                if s == "2H2":
                    factor=1e-3
                elif s == "T":
                    factor=1e-6
                else:
                    factor=1
                try:
                    line = ax2.plot(plot_slabs[islab].NH, plot_slabs[islab].GetAbd(s, nH=plot_nH)*factor, label=slabel, 
                               color=colors2[i], linestyle=linestyles[islab])
                except KeyError:
                    pass
                if linestyles[islab] == "-.":
                    line[0].set_dashes([10, 4, 2.2, 4])
                if linestyles[islab] == "--":
                    line[0].set_dashes([10, 5, 10, 5])
                if linestyles[islab] == ":":
                    line[0].set_dashes([2.5, 4, 2.5, 4])

        ax.set_xlim([5e18/Zd, 1e22/Zd]);
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend(fontsize=16.5, loc=3)
        ax.set_xlabel("$N/\mathrm{cm^{-2}}$", fontsize=20)
        ax.set_ylabel("$x_s = n_s/n$", fontsize=20)
        ax.set_ylim([1e-12*Zd,2e-3])
        ax.xaxis.set_ticks_position('bottom')
        ax3 = ax.twiny()
        ax3.set_xlim(np.array(ax.get_xlim())/1.87e21*Zd)
        ax3.set_xscale(ax.get_xscale())
        ax3.set_xlabel("$A_V$", fontsize=20)
        #ax2
        ax2.set_xlim(ax.get_xlim())
        ax2.set_ylim(ax.get_ylim())
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.legend(fontsize=16.5, loc=3)
        ax2.set_xlabel("$N/\mathrm{cm^{-2}}$", fontsize=20)
        #ax2.set_ylabel("$x_s = n_s/n$", fontsize=25)
        ax4 = ax2.twiny()
        ax4.set_xlim(np.array(ax.get_xlim())/1.87e21*Zd)
        ax4.set_xscale(ax.get_xscale())
        ax4.set_xlabel("$A_V$", fontsize=20)
        fig.tight_layout()
        fig.suptitle("$n_H={:d}$".format(plot_nH) + "$cm^{-3}$")
        if savefig:
            fig.savefig(fig_dir + figname + ".pdf")

#all Ni in appendix
def Plot_Ni(plot_slabs, plot_nH_arr, Zd, figname_end="", savefig=True):
    linestyles = ["-", "--", ":"]
    fig_dir = "/Users/munangong/chemistry_Athena/draft_chemistry/"

    if "OH" in figname_end:
        plot_spec_list=["CO", "C", "C+", "HCO+", "OHx", "CHx", "O", "O+", "OH+", "H2O+"]
    else:
        plot_spec_list=["CO", "C", "C+", "HCO+", "OHx", "CHx", "O", "O+"]
    plot_spec_label=["$\mathrm{CO}$", "$\mathrm{C}$", "$\mathrm{C^+}$", "$\mathrm{HCO^+}$", 
                     "$\mathrm{OH_x}$", "$\mathrm{CH_x}$", "$\mathrm{O}$", 
                     "$\mathrm{O^+}$", "$\mathrm{OH^+}$", "$\mathrm{H_2O^+}$"]
    colors = ["m", "r", "y", "g", "b", "cyan", "k", "gray", "BurlyWood", "DarkCyan"]

    plot_spec_list2=["2H2", "H+", "H2+", "H3+", "He+", "Si+", "e"]

    plot_spec_label2=[r"$2\times 10^{-3}\mathrm{H_2}$", "$\mathrm{H^+}$", "$\mathrm{H_2^+}$",
                      "$\mathrm{H_3^+}$", "$\mathrm{He^+}$",
                      "$\mathrm{Si^+}$", "$\mathrm{e}$"]
    colors2 = ["k", "y", "g", "b", "r", "cyan", "m"]

    for plot_nH in plot_nH_arr:
        plot_nH = int(plot_nH)
        figname = "Ni_all_nH"+str(plot_nH) + "Zd_" + str(int(Zd)) + figname_end
        fig=plt.figure(figsize=[20,8])
        ax = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        for islab in np.arange(len(plot_slabs)):
            for s,i in zip(plot_spec_list, np.arange(len(plot_spec_list))):
                if islab == 0:
                    slabel = plot_spec_label[i]
                else:
                    slabel = None
                try:
                    line = ax.plot(plot_slabs[islab].NH, 
                                   GetNi(plot_slabs[islab].NH, plot_slabs[islab].GetAbd(s, nH=plot_nH)), label=slabel, 
                                   color=colors[i], linestyle=linestyles[islab])
                    if linestyles[islab] == "-.":
                        line[0].set_dashes([10, 4, 2.2, 4])
                    if linestyles[islab] == "--":
                        line[0].set_dashes([10, 5, 10, 5])
                    if linestyles[islab] == ":":
                        line[0].set_dashes([2.2, 4, 2.2, 4])
                except KeyError:
                    pass
                if s == "CHx" or s == "OHx":
                    try:
                        line = ax.plot(plot_slabs[islab].NH, 
                                       GetNi(plot_slabs[islab].NH, plot_slabs[islab].GetAbd(s[:-1], nH=plot_nH)), 
                                       color=colors[i], linestyle="-.")
                        line[0].set_dashes([10, 4, 2.2, 4])
                    except KeyError:
                        pass
            for s,i in zip(plot_spec_list2, np.arange(len(plot_spec_list2))):
                if islab == 0:
                    slabel = plot_spec_label2[i]
                else:
                    slabel = None
                if s == "2H2":
                    factor=1e-3
                else:
                    factor=1
                try:
                    line = ax2.plot(plot_slabs[islab].NH, 
                                   GetNi(plot_slabs[islab].NH, plot_slabs[islab].GetAbd(s, nH=plot_nH))*factor, label=slabel, 
                                   color=colors2[i], linestyle=linestyles[islab])
                    if linestyles[islab] == "-.":
                        line[0].set_dashes([10, 4, 2.2, 4])
                    if linestyles[islab] == "--":
                        line[0].set_dashes([10, 5, 10, 5])
                    if linestyles[islab] == ":":
                        line[0].set_dashes([2.2, 4, 2.2, 4])
                except KeyError:
                    pass

        ax.set_xlim([5e19/Zd, 1e22/Zd]);
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend(fontsize=16.5, loc=3)
        ax.set_xlabel("$N/\mathrm{cm^{-2}}$", fontsize=20)
        ax.set_ylabel("$N_i/\mathrm{cm^{-2}}$", fontsize=20)
        ax.xaxis.set_ticks_position('bottom')
        ax.set_ylim([1e10*Zd, 1e19])
        ax3 = ax.twiny()
        ax3.set_xlim(np.array(ax.get_xlim())/1.87e21*Zd)
        ax3.set_xscale(ax.get_xscale())
        ax3.set_xlabel("$A_V$", fontsize=20)
        #ax2
        ax2.set_xlim(ax.get_xlim())
        ax2.set_ylim(ax.get_ylim())
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.legend(fontsize=16.5, loc=3)
        ax2.set_xlabel("$N/\mathrm{cm^{-2}}$", fontsize=20)
        ax4 = ax2.twiny()
        ax4.set_xlim(np.array(ax.get_xlim())/1.87e21*Zd)
        ax4.set_xscale(ax.get_xscale())
        ax4.set_xlabel("$A_V$", fontsize=20)
        fig.tight_layout()
        if savefig:
            fig.savefig(fig_dir + figname + ".pdf")

def plot_rates_nH(slab, nH, rate_spec, abd_spec=None, fig_size_factor=1.):
    if abd_spec is None:
        abd_spec = rate_spec
    rates_list_plot = slab.rates_list
    crates, drates = GetReaction(rates_list_plot, rate_spec)
    fig=plt.figure(figsize=[8*fig_size_factor,5.5*fig_size_factor])
    ax = fig.add_subplot(111)
    linestyles = ["-", "--", "-.", ":"]
    color_c = ["r", "orange", "yellow", "m"]
    color_d = ["b", "cyan", "g", "gray"]

    li=0
    ci = 0
    for r in crates:
        try:
            ax.plot(slab.NH, slab.GetAbd(r, nH=nH), 
                    label=r, color=color_c[ci], linestyle=linestyles[li],
                    linewidth=3)
        except KeyError:
            pass
        li += 1
        if li == 4:
            li=0
            ci += 1

    li=0
    ci = 0
    for r in drates:
        try:
            ax.plot(slab.NH, slab.GetAbd(r, nH=nH), 
                    label=r, color=color_d[ci], linestyle=linestyles[li],
                    linewidth=1.5)
        except KeyError:
            pass
        li += 1
        if li == 4:
            li=0
            ci += 1

    ax.set_xlim([1e19, 1e22]);
    ax2 = ax.twinx()
    try:
        ax2.plot(slab.NH, slab.GetAbd(abd_spec, nH=nH), color="k", linestyle="-")
    except KeyError:
        pass
    ax.set_xscale("log")
    ax2.set_xscale("log")
    ax.set_yscale("log")
    ax2.set_yscale("log")
    ax.set_ylim([1e-25, 1e-10])
    ax2.set_xlim([1e19/slab.Zd, 1e22/slab.Zd]);
    ax.legend(fontsize=10)
    ax.set_xlabel("$N_H/cm^{-2}$")
    ax.set_ylabel("rates / $s^{-1}$")
    ax2.set_ylabel("x({})".format(abd_spec))

    ax3 = ax2.twiny()
    ax3.set_xlim(np.array(ax2.get_xlim())/1.87e21*slab.Zd)
    ax3.set_xscale(ax2.get_xscale())
    ax3.set_xlabel("Av")
    ax.xaxis.set_ticks_position('bottom')

    ax.set_title(rate_spec+" reactions, nH={:.1f}\n\n\n".format(nH), fontsize=10)
    return

def plot_rates_NH(slab, NH, rate_spec, abd_spec=None, fig_size_factor=1.):
    if abd_spec is None:
        abd_spec = rate_spec
    rates_list_plot = slab.rates_list
    crates, drates = GetReaction(rates_list_plot, rate_spec)
    fig=plt.figure(figsize=[8*fig_size_factor,6*fig_size_factor])
    ax = fig.add_subplot(111)
    linestyles = ["-", "--", "-.", ":"]
    color_c = ["r", "orange", "yellow", "m"]
    color_d = ["b", "cyan", "g", "gray"]

    li=0
    ci = 0
    for r in crates:
        try:
            ax.plot(slab.nH, slab.GetAbd(r, NH=NH), 
                    label=r, color=color_c[ci], linestyle=linestyles[li],
                    linewidth=3)
        except KeyError:
            pass
        li += 1
        if li == 4:
            li=0
            ci += 1

    li=0
    ci = 0
    for r in drates:
        try:
            ax.plot(slab.nH, slab.GetAbd(r, NH=NH), 
                    label=r, color=color_d[ci], linestyle=linestyles[li],
                    linewidth=1.5)
        except KeyError:
            pass
        li += 1
        if li == 4:
            li=0
            ci += 1

    ax.set_xlim([1e19, 1e22]);
    ax2 = ax.twinx()
    try:
        ax2.plot(slab.nH, slab.GetAbd(abd_spec, NH=NH), color="k", linestyle="-")
    except KeyError:
        pass
    ax.set_xscale("log")
    ax2.set_xscale("log")
    ax.set_yscale("log")
    ax2.set_yscale("log")
    ax.set_ylim([1e-25, 1e-10])
    ax2.set_xlim([slab.nH[0], slab.nH[-1]]);
    ax.legend(fontsize=10)
    ax.set_xlabel("$n_H/cm^{-3}$")
    ax.set_ylabel("rates / $s^{-1}$")
    ax2.set_ylabel("x({})".format(abd_spec))
    ax.set_title(rate_spec+" reactions, NH={:.1e}\n".format(NH), fontsize=10)
    return ax

def plot_thermo_nH(ax, plot_slab, plot_nH,linestyle="-", label="", plot_unimportant = False, plot_legend=True):
    #important heating processes:
    ax.plot(plot_slab.NH, plot_slab.GetAbd("GPE", nH=plot_nH), color="orange", linestyle=linestyle, label="GPE, "+label)
    ax.plot(plot_slab.NH, plot_slab.GetAbd("GCR", nH=plot_nH), color="m", linestyle=linestyle, label="GCR")
    ax.plot(plot_slab.NH, plot_slab.GetAbd("Gtotal", nH=plot_nH), color="r", linestyle=linestyle, label="Gtotal", 
            linewidth=4)
    #important cooling processes:
    ax.plot(plot_slab.NH, plot_slab.GetAbd("LLya", nH=plot_nH), color="g",linestyle=linestyle, label="LLya")
    ax.plot(plot_slab.NH, plot_slab.GetAbd("LRec", nH=plot_nH), color="navy", linestyle=linestyle, label="LRec")
    ax.plot(plot_slab.NH, plot_slab.GetAbd("LOI", nH=plot_nH), color="lightblue", linestyle=linestyle, label="LOI")  
    ax.plot(plot_slab.NH, plot_slab.GetAbd("LCII", nH=plot_nH), color="cyan", linestyle=linestyle, label="LCII")
    ax.plot(plot_slab.NH, plot_slab.GetAbd("LCI", nH=plot_nH), color="olive", linestyle=linestyle, label="LCI")
    ax.plot(plot_slab.NH, plot_slab.GetAbd("LCOR", nH=plot_nH), color="lawngreen", linestyle=linestyle, label="LCOR")
    ax.plot(plot_slab.NH, plot_slab.GetAbd("Ltotal", nH=plot_nH), color="b", linestyle=linestyle, label="Ltotal", 
            linewidth=2.5)
    #unimportant processes:
    if plot_unimportant:
        ax.plot(plot_slab.NH, plot_slab.GetAbd("GH2diss", nH=plot_nH), color="y", linestyle="-", label="GH2diss")
        ax.plot(plot_slab.NH, plot_slab.GetAbd("GH2pump", nH=plot_nH), color="chocolate", linestyle="-", label="GH2pump")
        ax.plot(plot_slab.NH, plot_slab.GetAbd("GH2gr", nH=plot_nH), color="pink", linestyle="-", label="GH2gr")
        ax.plot(plot_slab.NH, plot_slab.GetAbd("LH2", nH=plot_nH), color="yellowgreen", linestyle="-", label="LH2")
        ax.plot(plot_slab.NH, plot_slab.GetAbd("LDust", nH=plot_nH), color="skyblue", linestyle="-", label="LDust")
        ax.plot(plot_slab.NH, plot_slab.GetAbd("LH2diss", nH=plot_nH), color="turquoise", linestyle="-", label="LH2diss")
        ax.plot(plot_slab.NH, plot_slab.GetAbd("LHIion", nH=plot_nH), color="indigo", linestyle="-", label="LHIion")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim([1e-30, 1e-23])
    ax.set_xlim([plot_slab.NH[0], plot_slab.NH[-1]])
    ax2 = ax.twinx()
    ax2.plot(plot_slab.NH, plot_slab.GetAbd("T", nH=plot_nH)/100, color="k",
            linestyle=linestyle, label="T/100", linewidth=4)
    ax2.plot(plot_slab.NH, plot_slab.GetAbd("e", nH=plot_nH), color="k", linestyle=linestyle, label="e",
            linewidth=2)
    ax2.plot(plot_slab.NH, plot_slab.GetAbd("H", nH=plot_nH), color="silver", linestyle=linestyle, label="H", 
            linewidth=3)
    ax2.plot(plot_slab.NH, plot_slab.GetAbd("H2", nH=plot_nH)*2, color="brown", linestyle=linestyle, label="2H2", 
            linewidth=1.5)
    ax2.plot(plot_slab.NH, plot_slab.GetAbd("C2Ctot", nH=plot_nH), color="gray", linestyle=linestyle,
            label="C/Ctot",linewidth=2)
    ax2.plot(plot_slab.NH, plot_slab.GetAbd("CO2Ctot", nH=plot_nH), color="brown", linestyle=linestyle,
            label="CO/Ctot", linewidth=3)
   
    ax2.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_ylim([1e-4, 2])
    ax2.set_xlim([plot_slab.NH[0], plot_slab.NH[-1]])
    if plot_legend:
        ax.legend(fontsize=10, loc=3)
        ax2.legend(fontsize=10, loc=2)
    ax.set_xlabel("$N_H/cm^{-2}$")
    ax2.set_ylabel("$T2, 2x_{H_2}, x_{CO}/x_{Ctot}, x_{C^+}/x_{Ctot}$")
    ax.set_ylabel("$\Lambda(\Gamma) / \mathrm{erg^{}s^{-1}H^{-1}}$")
    ax.set_title("Heating and cooling, nH={:.1f}".format(plot_nH), fontsize=10)
    return
