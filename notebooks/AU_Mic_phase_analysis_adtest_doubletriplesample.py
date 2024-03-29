from funcs.notebook import *
from funcs.phaseanalysis import (get_observed_phases,
                                 get_flare_phases,
                                 get_cumulative_distributions,)

from funcs.ad import sample_AD_for_custom_distribution, get_pvalue_from_AD_statistic
from scipy.stats import ks_1samp
from scipy import interpolate

from matplotlib.lines import Line2D

import datetime

import sys

ROTPER =  4.862 # martioli
ORBPER =  8.463 # exoplanet.eu 


def analyse_phase_distribution(subsample, sector, data, tstamp, mode, rotper=ROTPER, phaseshifts=[0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,
                       0.05,.15,.25,.35,.45,.55,.65,.75,.85,.95]):
    """Wrapper for the analysis of flare phases. 
    Generates the KS test value results.
    
    Parameters:
    -----------
    subsample
    sector
    data
    tstamp
    """
    # Select either one or both Sectors, and fetch their cadences
    aumic_, get_secs_cadences = data[sector]

    # Fetch the orbital phases for AU Mic b
    if mode != "Orbit":
        aumic_ = get_flare_phases(aumic_, "Rotation", rotper=rotper).reset_index()
    elif mode == "Orbit":
        aumic_ = get_flare_phases(aumic_, mode).reset_index()
    
  
    if subsample == "random half":
        
        indices = np.random.choice(aumic_.index.values, size=aumic_.shape[0]//2, replace=False)
        aumic = aumic_.iloc[indices,]
#         print(indices, indices.shape[0], aumic_.shape[0], aumic.shape[0])
        
    if subsample == "low energy half":
        aumic__ = aumic_[aumic_.ed_rec < 0.93] # ED = 0.93 s splits the Sector 1 sample in two halves
        
    if subsample == "high energy half":
        aumic__ = aumic_[aumic_.ed_rec > 0.93] # ED = 0.93 s splits the Sector 1 sample in two halves
        
    if subsample == "total":
        aumic__ = aumic_
        
    # Introduce artificial phase shift
    for phaseshift in phaseshifts:#
        print(f"Shift phase by {phaseshift}.")
        
        aumic = aumic__.copy()
        
        # shift phases    
        #print(aumic.columns)
        aumic.phases = (aumic.phases + phaseshift) % 1

        # Finally sort the flares by their phases in ascending order
        aumic = aumic.sort_values(by="phases", ascending=True)

        # define the "bins" as phase differences between the flare phases
        p = np.sort(aumic.phases.values)

	# SUPPLEMENT WITH double 

        # Get the observing times
        if mode !="Orbit":
            aumicphases, binmids = get_observed_phases("Rotation", p, get_secs_cadences,
                                                        rotper=rotper, phaseshift=phaseshift, test="KS")
        elif mode == "Orbit":
            aumicphases, binmids = get_observed_phases(mode, p, get_secs_cadences,
                                                      phaseshift=phaseshift, test="KS")
    
         # shift phases here too. NO THEY WERE SHIFTED ALREADY, see line above
#         aumicphases =  aumicphases + phaseshift

        # Get the (cumulative) flare phase distributions 
        n_i, n_exp, cum_n_exp, cum_n_i = get_cumulative_distributions(aumic, aumicphases, get_secs_cadences)

        # Define a cumulative cdf to pass to the K-S test

        # add the (0,0) and (1,1) points to the cdf
        cphases = np.insert(aumic.phases.values, 0, 0)
        cphases = np.insert(cphases, -1, 1)
        vals = np.insert(cum_n_exp, 0, 0)
        vals = np.insert(vals, -1, 1)

        # Interpolate!
        f = interpolate.interp1d(cphases, vals, fill_value="extrapolate")

        # plot diagnostic with the interpolated cdf for the AD test
        plt.figure(figsize=(8, 7))
        x = np.linspace(0, 1, 500)
        y = f(x)
        plt.scatter(aumic.phases, cum_n_exp, s=20, c="k", marker="x", label="expected distribution")
        plt.scatter(aumic.phases, cum_n_i, s=30, c="r", label="measured distribution")
        plt.plot(x, y, label="interpolation", zorder=-10, c="grey");
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.legend(loc=2, frameon=False)
        plt.xlabel("phase")
        plt.ylabel("cumul. fraction of observed flares")
        tst = tstamp.replace("-","_")
        plt.savefig(f"../results/plots/{tst}_AUMic_AD_Test_cumdist_{subsample}_{sector}_{mode}_shift{phaseshift}_triple.png", dpi=300)

        # double the sample size
        supp = np.random.normal((p[:-1] + p[1:]) / 2., np.abs(np.diff(p)) / 2., len(p) - 1)
        supp = np.insert(supp, -1, np.random.normal((p[-1] + p[0] + 1) / 2., (p[0] + 1 - p[-1]) / 2) % 1.)
        newp = np.sort(np.append(p, supp))
        
        # triple the sample size
#        newx1 = np.random.normal((p[:-1] * 2 + p[1:])/3.,np.abs(np.diff(p)) / 3.,len(p)-1)
#        newx2 = np.random.normal((p[:-1] + p[1:] * 2)/3.,np.abs(np.diff(p)) / 3.,len(p)-1)
#        newp = np.sort(np.concatenate([p, newx1, newx2]))
        print(newp.shape[0], p.shape[0], supp.shape[0])
        
        # Finally, the A-D tes
        N = 200000
        A2  = sample_AD_for_custom_distribution(f, newp.shape[0], N)
        #A2old = A2old[np.isfinite(A2old)]
        #A2extra = np.loadtxt(f"../results/{tst}_AUMic_AD_Test_cumdist_{subsample}_{sector}_{mode}_shift{phaseshift}_50000_triple_A2.txt")

        #A2 = np.concatenate([A2old, A2extra])
       # print(A2.shape)
        with open(f"../results/{tst}_AUMic_AD_Test_cumdist_{subsample}_{sector}_{mode}_shift{phaseshift}_{N}_double_A2.txt", "w") as file:
            np.savetxt(file, A2) 

        pval, atest = get_pvalue_from_AD_statistic(newp, f, A2)
        print(pval, atest)

        # plot diagnostic with the A2 distribution
        plt.figure(figsize=(8, 7))
        
        plt.hist(A2, bins=np.linspace(np.min(A2),np.max(A2),N//1000))
        plt.axvline(atest,c="k")    
        plt.xlim(np.min(A2),np.max(A2))
        
        plt.xlabel("A2")
        tst = tstamp.replace("-","_")
        plt.savefig(f"../results/plots/{tst}_AUMic_AD_Test_A2dist_{subsample}_{sector}_{mode}_shift{phaseshift}_double.png", dpi=300)


        

        # write out results
        with open("../results/adtests.csv", "a") as f:
            #tstamp	period	sector	subsample	AD	p	nsteps_mcmc	nflares	totobs_days	shift
            stri = (f"{tstamp},{mode},{sector},{subsample},{atest},"
                    f"{pval},{newp.shape[0]},{A2.shape[0]},"
                    f"{aumicphases.sum().sum()/60./24.},{phaseshift}\n")
            f.write(stri)
            
        # free up space for next iteration, just being nice to your cache
        del aumic
        
def paper_figure_kstest(tstamp):
    
    # Get the data
    kss_ = pd.read_csv("../results/adtests.csv").drop_duplicates()
    kss_ = kss_[kss_.tstamp == tstamp]
    
    
    # Set up figure
    plt.figure(figsize=(14,7))

    # Rotation to the left
    plt.subplot(121)

    kss = kss_[kss_.period == "Rotation"]

    colors = {"Both Sectors" : ("k","o"),
              "Sector 27" : ("r","d"),
              "Sector 1" : ("c","X")}

    handles = []
    for label, g in kss.groupby("sector"):    

        c, m = colors[label]
        for ll, h in g.groupby("phaseshift"):
            
            plt.scatter(h.subsample, h["p"], marker=fr"${ll}$", s=240, c=c, alpha=1)
        handles.append((Line2D([0], [0], marker=m, color='w', 
                               label=label, markerfacecolor=c, 
                               markersize=15)))

    # make sigma thresholds    
    plt.axhline(1 - .342*2, c="k", linestyle="dashed")
    plt.axhline(1 - .342*2 - .136*2, c="k", linestyle="dotted")
    plt.axhline(1 - .342*2 - .136*2 - .021*2, c="k", linestyle="dashdot")


    plt.title("Rotational modulation")

    plt.ylim(1e-3,1)
    plt.yscale("log")
    plt.ylabel(r"$p$ value")

    plt.legend(handles=handles, frameon=False, loc=(0.07,0.1));

    # Orbit to the right
    plt.subplot(122)

    kss = kss_[kss_.period == "Orbit"]

    handles = []
    for label, g in kss.groupby("sector"):
        c, m = colors[label]
        if label != "Both Sectors":
            zorder = -10
        else:
            zorder=1
        for ll, h in g.groupby("phaseshift"):
            plt.scatter(h.subsample, h["p"],  marker=fr"${ll}$", s=240,
                        c=c, alpha=1, zorder=zorder)
        handles.append((Line2D([0], [0], marker=m, color='w', 
                               label=label, markerfacecolor=c, 
                               markersize=15)))

    plt.legend(handles=handles, frameon=False, loc=(0.3,0.1))

    # make sigma thresholds    
    plt.axhline(1 - .342*2, c="k", linestyle="dashed")
    plt.axhline(1 - .342*2 - .136*2, c="k", linestyle="dotted")
    plt.axhline(1 - .342*2 - .136*2 - .021*2, c="k", linestyle="dashdot")

    plt.title("Orbital modulation")
    plt.ylim(1e-3,1)
    plt.yscale("log")
    plt.text(2.2,.3,r"1 $\sigma$")
    plt.text(2.2,.045,r"2 $\sigma$")
    plt.text(2.2,.002,r"3 $\sigma$")


    plt.tight_layout()
    tst = tstamp.replace("-","_")
    plt.savefig(f"../results/plots/{tst}_AUMic_KStests_meta.png", dpi=300)
    plt.savefig(f"/home/ekaterina/Documents/002_writing/aumic-flaring-spi-draft/figures/{tst}_AUMic_ADtests_meta.png",
                dpi=300)

def paper_figure_kstest_mode(tstamp, mode):
    
    # Get the data
    kss_ = pd.read_csv("../results/adtests.csv").drop_duplicates()
    kss_ = kss_[kss_.tstamp == tstamp]
    
    # Set up figure
    plt.figure(figsize=(7,7))


    kss = kss_[kss_.period == mode]

    colors = {"Both Sectors" : ("k","o"),
              "Sector 27" : ("r","d"),
              "Sector 1" : ("c","X")}

    handles = []
    for label, g in kss.groupby("sector"):    
        c, m = colors[label]
        for ll, h in g.groupby("shift"):
            
            plt.scatter(h.subsample, h["p"], marker=fr"${ll}$", s=240, c=c, alpha=1)
        handles.append((Line2D([0], [0], marker=m, color='w', 
                               label=label, markerfacecolor=c, 
                               markersize=15)))
#         plt.scatter(g.subsample, g["p"], marker=m, s=140, c=c)
#         handles.append((Line2D([0], [0], marker=m, color='w', 
#                                label=label, markerfacecolor=c, 
#                                markersize=15)))

    # make sigma thresholds    
    plt.axhline(1 - .342*2, c="k", linestyle="dashed")
    plt.axhline(1 - .342*2 - .136*2, c="k", linestyle="dotted")
    plt.axhline(1 - .342*2 - .136*2 - .021*2, c="k", linestyle="dashdot")


    plt.title(f"{mode} modulation")

    plt.ylim(1e-5,1)
#    plt.yscale("log")
    plt.ylabel(r"$p$ value")

    plt.legend(handles=handles, frameon=False, loc=(0.55,0.1));


    plt.tight_layout()
    tst = tstamp.replace("-","_")
    plt.savefig(f"../results/plots/{tst}_AUMic_ADtests_meta_{mode}.png", dpi=300)
#    plt.savefig(f"/home/ekaterina/Documents/002_writing/aumic-flaring-spi-draft/figures/{tst}_AUMic_ADtests_meta.png",
#                dpi=300)




if __name__ == "__main__":
    
    # SETUP
    sector = int(sys.argv[1])
    subs = int(sys.argv[2])

    # Load data and select the final vetted flares

    # Sector 1
    aumic1 = pd.read_csv("../results/2021_02_18_AUMic_flares_1.csv")
    aumic1 = aumic1[(aumic1.final==1) & (aumic1["real?"]==1)]

    # Sector 27
    aumic27 = pd.read_csv("../results/2021_02_11_AUMic_flares_27.csv")
    aumic27 = aumic27[(aumic27.final==1) & (aumic27["real?"]==1)]

    # Both Sectors
    aumic127 = pd.concat([aumic1, aumic27])

    # Define the sector number and observing cadence for each light curve and link to selection
    data = {"Both Sectors": [aumic127,[(1, 2), (27, 1/3)]], 
            "Sector 1": [aumic1, [(1, 2)]], 
            "Sector 27": [aumic27, [(27, 1/3)]]}

    # list the subsamples you want to analyse
    subsamples = ["total", "high energy half", "low energy half"]
    # get the keys to loop
    sectors = list(data.keys())

    # timestamp for unique rows in results table
    tstamp = "2021-12-20"#datetime.date.today().isoformat()

    
    # ANALYSIS
    
    # Loop over all configurations
#    mode = "Orbit"
#    for subsample in subsamples[:1]:
#        for sector in sectors[1:2]:
#            print(f"{mode}: Analyzing sumbsample {subsample}, Sector {sector}")
#            analyse_phase_distribution(subsample, sector, data, tstamp, mode)

  

    if sys.argv[3]=="b":
        mode = "Beat Period"
        per =  1. / ((1. / ROTPER) - (1. / ORBPER)) # martioli
    elif sys.argv[3]=="r":
        mode = "Rotation"
        per = ROTPER
    elif sys.argv[3]=="o":
        mode = "Orbit"
        per = ROTPER
#    for subsample in subsamples[2:]:
#        for sector in sectors[:]: #DO loop over both Sectors in Rotation mode
#            print(f"{mode}: Analyzing sumbsample {subsample}, Sector {sector}")
#            analyse_phase_distribution(subsample, sector, data, tstamp, mode, rotper=per)

    phaseshifts = [float(sys.argv[4])]

    print(f"phase shift {phaseshifts}")

    analyse_phase_distribution(subsamples[subs], sectors[sector], data, tstamp, mode, rotper=per, phaseshifts=phaseshifts)
            
 #   analyse_phase_distribution("total", "Both Sectors", data, tstamp, "Orbit", rotper=ROTPER)

    # Generate paper figure
#    paper_figure_kstest_mode(tstamp, "Beat Period")
#    paper_figure_kstest_mode(tstamp, "Rotation")
#    paper_figure_kstest_mode(tstamp, "Orbit")
