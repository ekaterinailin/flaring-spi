from funcs.notebook import *
from funcs.phaseanalysis import (get_observed_phases,
                                 get_flare_phases,
                                 get_cumulative_distributions,
                                )


from funcs.ad import (anderson_darling_statistic,
                      get_pvalue_from_AD_statistic,
                      sample_AD_for_custom_distribution,
                      )

#from scipy import interpolate

from matplotlib.lines import Line2D

import datetime
import time

ROTPER =  4.862 # martioli
ORBPER =  8.463 # exoplanet.eu


def analyse_phase_distribution_ad(subsample, sector, data, tstamp, mode, N, rotper=ROTPER, phaseshift=0):
    """Wrapper for the analysis of flare phases. 
    Generates the KS test value results.
    
    Parameters:
    -----------
    subsample
    sector
    data
    tstamp
    N  - MCMC steps
    """
    # Select either one or both Sectors, and fetch their cadences
    aumic_, get_secs_cadences = data[sector]

    # Fetch the orbital phases for AU Mic b
    if mode != "Orbit":
        aumic_ = get_flare_phases(aumic_, "Rotation", rotper=rotper).reset_index()
    else:
        aumic_ = get_flare_phases(aumic_, mode, rotper=rotper).reset_index()
  
    if subsample == "random half":
        indices = np.random.choice(aumic_.index.values, size=aumic_.shape[0]//2, replace=False)
        
        aumic = aumic_.iloc[indices,]
        print(indices, indices.shape[0], aumic_.shape[0], aumic.shape[0])
       
    elif subsample == "low energy half":
        aumic = aumic_[aumic_.ed_rec < 0.93] # ED = 0.93 s splits the Sector 1 sample in two halves
        
    elif subsample == "high energy half":
        aumic = aumic_[aumic_.ed_rec > 0.93] # ED = 0.93 s splits the Sector 1 sample in two halves
        
    elif subsample == "total":
        aumic = aumic_
        

    aumic.phases = (aumic.phases.values + phaseshift) % 1.
        
    # Finally sort the flares by their phases in ascending order
    aumic = aumic.sort_values(by="phases", ascending=True)

    # define the "bins" as phase differences between the flare phases
    p = aumic.phases.values

    # Get the observing times
    if mode!= "Orbit":
        aumicphases, ph = get_observed_phases("Rotation", np.linspace(0.01,.99,100), 
                                              get_secs_cadences, rotper=rotper, 
                                               phaseshift=phaseshift)
    else:
        aumicphases, ph = get_observed_phases(mode, np.linspace(0.01,.99,100), 
                                      get_secs_cadences, rotper=rotper, 
                                       phaseshift=phaseshift)
    # observing times weighted by number of flares per sector to account for detection threshold variability
    #cobs = aumicphases.sum(axis=1).values
    for qcs, cadence in get_secs_cadences:
#        print(aumicphases[qcs])
        aumicphases[qcs] = aumicphases[qcs] * aumic[aumic.qcs == qcs].shape[0] / aumic.shape[0]
#        print(aumicphases[qcs])
    cobs = aumicphases.sum(axis=1).values

    # sample from the near-uniform expected distribution and get the AD statistic

    A2 = sample_AD_for_custom_distribution(p, cobs, N, ph, l=phaseshift)
#    print(A2)
    
    # get p value and statistic value
    pvalue, adtestvalue = get_pvalue_from_AD_statistic(p, A2)
    
    # plot diagnostic with the interpolated callable cdf for the AD test
    plt.figure(figsize=(8, 7))
    plt.title(f"{subsample}, {sector}, {mode}")
    plt.hist(A2, bins=N//100, edgecolor="k", histtype="step", label="AD statistic distribution")
    plt.axvline(adtestvalue, label="measured AD")
    plt.legend(loc=1, frameon=False)
    plt.xlabel("AD")
    plt.ylabel("counts")
    tst = tstamp.replace("-","_")
    plt.savefig(f"../results/plots/{tst}_AUMic_AD_Test_{subsample}_{sector}_{mode}_{phaseshift}.png", dpi=300)
    
    # write out results
    with open("../results/adtests.csv", "a") as f:
        #tstamp	period	sector	subsample	D+	p+	D-	p-	D	p	nflares	totobs_days
        stri = (f"{tstamp},{mode},{sector},{subsample},"
                f"{adtestvalue},{pvalue},{N},{aumic.shape[0]},"
                f"{aumicphases.sum().sum()/60./24.},{phaseshift}\n")
        f.write(stri)
        
def paper_figure_adtest(tstamp):
    
    # Get the data
    kss_ = pd.read_csv("../results/adtests.csv").drop_duplicates()
    kss_ = kss_[kss_.tstamp == tstamp]
    kss_ = kss_[kss_.nsteps_mcmc == 10000]
    
    
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


    plt.title("Rotational modulation")

    plt.ylim(1e-3,1)
    plt.yscale("log")
    plt.ylabel(r"$p$ value")

    plt.legend(handles=handles, frameon=False, loc=(0.5,0.1));

    # Orbit to the right
    plt.subplot(122)

    kss = kss_[kss_.period == "Orbit"]


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

    plt.legend(handles=handles, frameon=False, loc=(0.6,0.1))

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
    plt.savefig(f"../results/plots/{tst}_AUMic_ADtests_meta.png", dpi=300)
    plt.savefig(f"/home/ekaterina/Documents/002_writing/aumic-flaring-spi-draft/figures/{tst}_AUMic_ADtests_meta.png",
                dpi=300)


def paper_figure_adtest_mode(tstamp, mode):
    
    # Get the data
    kss_ = pd.read_csv("../results/adtests.csv").drop_duplicates()
    kss_ = kss_[kss_.tstamp >= '2021-09-01']
    print(kss_.shape)
    kss_ = kss_[kss_.nsteps_mcmc == 10000]
    
    
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

    plt.ylim(1e-3,1)
    plt.yscale("log")
    plt.ylabel(r"$p$ value")

    plt.legend(handles=handles, frameon=False, loc=(0.5,0.1));


    plt.tight_layout()
    tst = tstamp.replace("-","_")
    plt.savefig(f"../results/plots/{tst}_AUMic_ADtests_meta_{mode}.png", dpi=300)
#    plt.savefig(f"/home/ekaterina/Documents/002_writing/aumic-flaring-spi-draft/figures/{tst}_AUMic_ADtests_meta.png",
#                dpi=300)



if __name__ == "__main__":
    
    # SETUP

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
    
    # list phaseshifts
#    shifts = [.0, .2, .4, .6, .8,]
    shifts = [.1, .3, .5, .7, .9,]

    # timestamp for unique rows in results table
    tstamp = datetime.date.today().isoformat()

    
    # ANALYSIS
    
    
    # Loop over all configurations
#    mode = "Orbit"
#    for subsample in subsamples[1:]:
#        for sector in sectors:
#            for shift in shifts:
#                print(f"{mode}: Analyzing sector {sector}, shift {shift}, subsample {subsample}")
#                analyse_phase_distribution_ad(subsample, sector, data, tstamp, mode, 10000, phaseshift=shift)

#     mode = "Rotation"
#     ROTPER =  4.862 # martioli
#     for subsample in subsamples[:1]:
#         for sector in sectors[:]: #DO loop over both Sectors in Rotation mode
#             for shift in shifts:
#                 print(f"{mode}: Analyzing sector {sector}, shift {shift} subsample {subsample}")
#                 analyse_phase_distribution_ad(subsample, sector, data, tstamp, mode, 10000, phaseshift=shift)


    # Generate paper figure
#    paper_figure_adtest(tstamp)



    mode = "Beat Period"
    per =  1. / ((1. / ROTPER) - (1. / ORBPER)) # martioli
    print(per)
    for subsample in subsamples[:1]:
        for sector in sectors[:1]: #DO loop over both Sectors in Rotation mode
            for shift in shifts:
                print(f"{mode}: Analyzing sector {sector}, shift {shift} subsample {subsample}")
                analyse_phase_distribution_ad(subsample, sector, data, tstamp, 
                                              mode, 10000, phaseshift=shift, rotper=per)
#                time.sleep(20)

#    mode = "Beat Period"
#    paper_figure_adtest_mode(tstamp, mode)

#    mode = "Rotation"
#    paper_figure_adtest_mode(tstamp, mode)

#    mode = "Orbit"
#    paper_figure_adtest_mode(tstamp, mode)

                
