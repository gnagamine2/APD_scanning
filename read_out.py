    # -*- coding: utf-8 -*-
"""
original file created by: rober
adapted for BSc Thesis by Julian Santen (from Sept. 2022)
"""

# Import Custom Function
from functionDefinitions.functions_read_out import *


def create_folder(filepath):
    if not (os.path.exists(filepath)): # if folder not found
        os.mkdir(filepath)
        print("folder created: ", filepath)
    else:
        #print("folder found: ", filepath)
        pass

#%% MAIN FUNCTION


if __name__ == "__main__":

    ########################################################
    # Data Import Here --> Need to change paths
    ########################################################

    # parameters (use forward slashes!)
    """
    folder = "/Users/JulianSanten/GIT/OMEL_GIT/Data/23-09-23_AP-3-108_TestingJulianScript_IntesityTracingLargeArea+Spec/"
    #namelist = ["AP-3-108_Conc10-6_try1x_260_y108_z97_range30_pts150_Particle_0_x273y_100"]
    #namelist = ["AP-3-108_Conc10-6_try1x_260_y108_z97_range30_pts150_Particle_1_x262y_99"]
    namelist = ["AP-3-108_10-6conc_1MHz_SingleDotTry"]
    file_extension_is_out = False
    is_Aniket = False
    """

    """
    folder = '/Users/JulianSanten/GIT/OMEL_GIT/read-tttr/Data/'
    namelist = ['Hydraharp_2022-06-02_PM_171'] #shows very nice g(2) characteristics
    file_extension_is_out = False
    is_Aniket = False
    """

    """
    folder = "/Users/JulianSanten/GIT/OMEL_GIT/Data/22-09-01_ScanningStageTTTRCodeTest3/"
    namelist = ["Test_AP_4_17_10-6_RepRate2,5x10^-6_x_274_y60_z94_range30_pts150_Particle_0_x281y_70"]
    
    #namelist = ["Test_AP_4_17_10-6_RepRate2,5x10^-6_x_274_y60_z94_range30_pts150_Particle_0_x281y_70", "Test_AP_4_17_10-6_RepRate2,5x10^-6_x_274_y60_z94_range30_pts150_Particle_0_x288y_49",
    #                "Test_AP_4_17_10-6_RepRate2,5x10^-6_x_274_y60_z94_range30_pts150_Particle_1_x279y_58", "Test_AP_4_17_10-6_RepRate2,5x10^-6_x_274_y60_z94_range30_pts150_Particle_1_x287y_72"]
    
    file_extension_is_out = True #Set if the file extension is .out or .ptu
    is_Aniket = False
    """

    #Anikets Dataset
    """
    folder = "/Users/JulianSanten/GIT/OMEL_GIT/Data/20200723_AP_109_oldsample_g2/"
    namelist = ["AP_109_old_P2_10MHz"] #great Decay!
    
    #namelist = ["AP_109_old_P1_10MHz_meas2", "AP_109_old_P1_10MHz_meas3_1ODless", "AP_109_old_P1_10MHz",
    #            "AP_109_old_P3_10MHz_differentND_ND1dropduringmeas", "AP_109_old_P2_10MHz",
    #            "AP_109_old_P4_10MHz_differentND4_ND1dropduringmeas", "AP_109_old_P5_10MHz_ND3", "AP_109_old_P6_10MHz_ND3"]
    file_extension_is_out = True
    is_Aniket = True #Aniket used a different naming scheme
    """

    """
    folder = "C:/Users/omel-acton/Desktop/Gabriel_Local_EMCCD/22-11-08/"
    namelist = ["AP_3_108_Run3_1e-7_x_150_y130_z95_range30_pts150_Particle_0_x163y_142"]
    file_extension_is_out = True
    is_Aniket = False #Aniket used a different naming scheme
    """

    """
    folder = "/Users/JulianSanten/GIT/OMEL_GIT/Data/22-11-07-dry-run-debugging/"
    # namelist = ["DryRun6_x_183_y161_z96_range30_pts150_Particle_0_x192y_162"]
    namelist = []
    file_extension_is_out = True
    is_Aniket = False  # Aniket used a different naming scheme
    """
    

    folder = ["C:/Users/gabri/Documents/GitHub/APD_scanning/Data/22-11-17-AP-3-120_MultipleSingleDot_Measurement/"]
    #folder = "T:/Gabriel/Data/22-11-17-AP-3-120_MultipleSingleDot_Measurement/"

    
    """
    # 0 Dots AP_4_24
    folder = "T:/Gabriel/Data/22-11-23/AP_4_24/FinalMeasurement/"
    """

    """
    # 50 Dots AP_3_108
    folder = "T:/Gabriel/Data/22-11-23/AP_3_108/"
    """
    # All files for BSc JSANTEN
    
    list_folders = ["C:/Users/gabri/Documents/GitHub/APD_scanning/Data/22-11-17-AP-3-120_MultipleSingleDot_Measurement/"]
    
    # list_folders = ["T:/Gabriel/Data/22-11-08_AP-3-108_MultipleDots_2Batches/",
    #                 "T:/Gabriel/Data/22-11-10-AP-3-120_MultipleSingleDot_Measurement/",
    #                 "T:/Gabriel/Data/22-11-11-AP-3-108_MultipleSingleDot_Measurement/",
    #                 "T:/Gabriel/Data/22-11-15-AP-4-17_MultipleSingleDot_Measurement/", 
    #                 "T:/Gabriel/Data/22-11-17-AP-3-120_MultipleSingleDot_Measurement/",
    #                 "T:/Gabriel/Data/22-11-18-AP-3-108_MultipleSingleDot_Measurement/",
    #                 "T:/Gabriel/Data/22-11-21-AP-4-17_and_AP-4-24try_MultipleSingleDot_Measurement-MeasurementFailMidleNight/AP-4-17/",
    #                 "T:/Gabriel/Data/22-11-23-AP-3-108_and_AP-4-24try_MultipleSingleDot_Measurement/AP_3_108/",
    #                 "T:/Gabriel/Data/22-11-23-AP-3-108_and_AP-4-24try_MultipleSingleDot_Measurement/AP_4_24/FinalMeasurement/"]

    # list_folders = ["C:/Users/gabri/Documents/GitHub/apd_scanning_imaging/Data/"]


    for folder in list_folders:

        namelist = []
        file_extension_is_out = True
        is_Aniket = False  # Aniket used a different naming scheme
    
        for file in os.listdir(folder + "measurements/"):
            if file.endswith(".out"):
                #print(os.path.join(folder, file[:-13]))
                namelist.append(file[:-13])
        """
        namelist = ["AP-3-120_C1x10-6_SingleDot_2,5x10+6Hz_405exct_250nW_ptu_x_77_y93_z124_range30_pts150_Particle_91_x74y_103",
                    "AP-3-120_C1x10-6_SingleDot_2,5x10+6Hz_405exct_250nW_ptu_x_77_y93_z124_range30_pts150_Particle_78_x77y_89"]
        
        namelist = ["AP-3-120_C1x10-6_SingleDot_2,5x10+6Hz_405exct_250nW_ptu_x_77_y93_z124_range30_pts150_Particle_78_x77y_89",
                    "AP-3-120_C1x10-6_SingleDot_2,5x10+6Hz_405exct_250nW_ptu_x_77_y93_z124_range30_pts150_Particle_24_x87y_103",
                    "AP-3-120_C1x10-6_SingleDot_2,5x10+6Hz_405exct_250nW_ptu_x_77_y93_z124_range30_pts150_Particle_92_x74y_106",
                    "AP-3-120_C1x10-6_SingleDot_2,5x10+6Hz_405exct_250nW_ptu_x_77_y93_z124_range30_pts150_Particle_75_x77y_103",
                    "AP-3-120_C1x10-6_SingleDot_2,5x10+6Hz_405exct_250nW_ptu_x_77_y93_z124_range30_pts150_Particle_131_x67y_78"]
        namelist = ["AP-3-120_C1x10-6_SingleDot_2,5x10+6Hz_405exct_250nW_ptu_x_77_y93_z124_range30_pts150_Particle_24_x87y_103"]
        #namelist = ["AP-3-120_C1x10-6_SingleDot_2,5x10+6Hz_405exct_250nW_ptu_x_77_y93_z124_range30_pts150_Particle_131_x67y_78"]
        namelist = ["AP-3-120_C1x10-6_SingleDot_2,5x10+6Hz_405exct_250nW_ptu_x_77_y93_z124_range30_pts150_Particle_78_x77y_89"]
        """
        #print("namelist: ", namelist)
        
        # Set Paramters
        binwidth = 0.1  # width of time bins to create intensity trace [s]
        nrbins = 1  # intensity bin to separate intensity trace in (nrbins) to separately evaluate different int. levels
        histbinmultiplier = 1  # (histbinmultiplier)*dtmicro = bin size for lifetime histogram
        dttau = 70e-9  # end time of lifetime fit [s]. Note: this is not the width but the end time!
    
        savebool = True
        singlebintaubool = False # IMPORTANT: Computation takes several minutes per blob
        g2bool = True
    
    
        ########################################################
        # LEAVE UNCHANGED
        ########################################################
    
        # Check if Necessary Folders are created (namely: /output/pdf)
        output_folder = os.path.join(folder, 'output')
        create_folder(output_folder)
    
        pdf_output_folder = os.path.join(output_folder, 'pdf')
        create_folder(pdf_output_folder)
    
        analysis_folder = os.path.join(folder, 'analysis/')
        create_folder(analysis_folder)
    
        g2_folder = os.path.join(analysis_folder, "g2/")
        create_folder(g2_folder)
        create_folder(os.path.join(g2_folder, "discarded/"))
        create_folder(os.path.join(g2_folder, "singular/"))
        create_folder(os.path.join(g2_folder, "non_singular/"))
    
        decay_folder = os.path.join(analysis_folder, "decay/")
        create_folder(decay_folder)
    
        int_trace_folder = os.path.join(analysis_folder, "int_trace/")
        create_folder(int_trace_folder)
    
        lifetime_folder = os.path.join(analysis_folder, "lifetime/")
        create_folder(lifetime_folder)
    
    
        # create arrays for export
        state_files = pd.DataFrame()
        no_cols_df = 0
    
    
        nrmeas = len(namelist)  # nr of data files
    
        for measnr in tqdm(range(nrmeas)):  # for each measurement file
            #%%  load data depending on the file extension
            print("")
            print("")
            print("file: ", namelist[measnr])
            if file_extension_is_out:
                if is_Aniket:
                    settingsfile = namelist[measnr] + '_settings'  # HHsettings file should have the same file name apart from this last identifier
                    HHsettings = load_obj(settingsfile, folder)
                    data = ImportT3(folder + namelist[measnr] + '.out', HHsettings)  # import data from .out file
    
                    [dtmicro, dtmacro, microtimes0, macrotimes0, microtimes1, macrotimes1, nrphotons0, nrphotons1, overflows,
                     microtimesP, macrotimesP, nrP, microtimesL, macrotimesL, nrL] = importT3marked(
                        folder + namelist[measnr] + '.out', HHsettings)  # import data from .out file
                else:
                    #Load HHSettings
                    settingsfile = namelist[measnr] + '_HHsettings'  # HHsettings file should have the same file name apart from this last identifier
                    HHsettings = load_obj("measurements/" + settingsfile, folder)
                    data = ImportT3(folder + "measurements/" + namelist[measnr] + '_tttrmode.out', HHsettings)  # import data from .out file
                    [dtmicro, dtmacro, microtimes0, macrotimes0, microtimes1, macrotimes1, nrphotons0, nrphotons1, overflows,
                     microtimesP, macrotimesP, nrP, microtimesL, macrotimesL, nrL] = importT3marked(
                        folder + "measurements/" + namelist[measnr] + '_tttrmode.out', HHsettings)  # import data from .out file

                no_cols_df = state_files.shape[0] # to add the data to the previous row
                state_files.at[no_cols_df, "name"] = namelist[measnr]
                state_files.at[no_cols_df, "nrphotons0"] = nrphotons0
                state_files.at[no_cols_df, "nrphotons1"] = nrphotons1
            else:
                data = ImportT3(folder + namelist[measnr] + '.ptu', 0) # import data from .ptu file
    
            # %%
            # macrotime values can store up to 2^10 = 1024 before overflow. overflows*1024*dtmacro gives experiment time [s]
            Texp = round(data[8] * data[1] * 1024,0)
            print('averaged cps on det0 and det1:', np.array(data[6:8]) / Texp)
            print('experimental time in s:', Texp)
    
            # %% init some arrays and values
            nrtbins = int(Texp / binwidth)  # nr of bins in intensity trace
            Ilimslist = np.ndarray(shape=(nrmeas, nrbins), dtype=float)
            tauavelist = np.full((nrmeas, nrbins), 0, dtype=float)  # init list to save lifetimes of each intensity bin
            Alist = np.full((nrmeas, nrbins), 0, dtype=float)  # init list to save maximum histogram amplitudes (at t=t0)
            taubinlist = np.zeros(nrtbins - 1, dtype=float)  # init list to save lifetimes of each time bin (binwidth)
            intbinlist = np.zeros(nrtbins - 1, dtype=float)  # init list to save intensity per bin
    
            # compensate detector offset
            [microtimes0, microtimes1, times0, times1, dtmicro, dtmacro, decaytlist, decayylist, fig] = ShiftPulsedData(data[2],
                                                                                                                   data[4],
                                                                                                                   data[3],
                                                                                                                   data[5],
                                                                                                                   data[0],
                                                                                                                   data[1])
            #plt.title(namelist[measnr] +  " Particle No: " + str(measnr))
            fig.tight_layout()
            #plt.show()
            if savebool == True:
                np.save(folder + 'output/' + namelist[measnr] + '_decay_tlist', decaytlist)
                np.save(folder + 'output/' + namelist[measnr] + '_decay_ylist', decayylist)
                try:
                    fig.savefig(decay_folder + namelist[measnr] + '_decay_tlist.pdf')
                    fig.savefig(decay_folder + namelist[measnr] + '_decay_tlist.png', dpi=300, transparent=False)
                except OSError as e:
                    print(e)
                    warnings.warn("Create the folder")
                    
            plt.close(fig)
    
            tmax = decaytlist[decayylist[:-1].argmax()]  # find histogram time with max photon count [ns]
    
            # get number of first photon per bin for both detectors (APD1 = 0, APD2 = 1)
            limits0 = HistPhotons(times0 * dtmicro, binwidth, Texp)  # gives index of first photon of det0 in each bin
            limits1 = HistPhotons(times1 * dtmicro, binwidth, Texp)
    
            # make an intensity trace and find Imax
            inttrace, fig, xIntTrace, yIntTrace = MakeIntTrace(limits0, limits1, binwidth, Texp)
            sns.set_palette("bright")
            #plt.suptitle(namelist[measnr])
            fig.tight_layout()
            #fig.show()
            if savebool == True:
                try:
                    fig.savefig(int_trace_folder + namelist[measnr] + '_IntTrace.pdf')
                    fig.savefig(int_trace_folder + namelist[measnr] + '_IntTrace.png', dpi=300, transparent=False)
                    np.savetxt(int_trace_folder + namelist[measnr] + '_csvIntTrace.csv', np.transpose([xIntTrace, yIntTrace]), delimiter = ",", header="x,y", comments="")
                except OSError as e:
                    print(e)
                    warnings.warn("Create the folder")
            plt.close(fig)
    
    
            Imax = np.max(inttrace)  # search overall maximum intensity per intensity bin
            Ilims = np.arange(0, Imax * 0.95 + 1, Imax * 0.95 / nrbins,
                              dtype=float)  # separate intensities into (nrbins) intensity bins
            Ilimslist[measnr, :] = np.array([Ilims[binnr] * 0.5 + Ilims[binnr + 1] * 0.5 for binnr in
                                             range(nrbins)])  # get centered intensity value per intensity bin
    
            # get G2
            if g2bool == True:
                fig = plt.figure(figsize=(10,5))
                [g2tlist, g2ylist, _, _, state, g2_tau_zero] = MakeG2(times0, times1, dtmicro, dtmacro, microtimes0, microtimes1,
                                                         g2restime=8e-9, nrbins=1000)
    
                fig.tight_layout()
                #plt.show()
    
                # store the results in a dataframe (for future use in the spectra script)
                state_files.at[no_cols_df, "state"] = state
                state_files.at[no_cols_df, "g2_tau_zero"] = g2_tau_zero
    
                if savebool == True:
                    np.save(folder + 'output/' + namelist[measnr] + '_G2_tlist', g2tlist)
                    np.save(folder + 'output/' + namelist[measnr] + '_G2_ylist', g2ylist)
                    try:
                        fig.savefig(g2_folder + state + "/" + namelist[measnr] + '_G2_ylist.pdf')
                        fig.savefig(g2_folder + state + "/" + namelist[measnr] + '_G2_ylist.png', dpi=300, transparent=False)
                    except OSError as e:
                        print(e)
                        warnings.warn("Create the g2 folder")
                plt.close(fig)
    
            for binnr in range(nrbins):  # for each intensity bin
                # select photons in bins with intensity within Ilims(binnr) and Ilims(binnr+1)
                [binmicrotimes0, bintimes0, binmicrotimes1, bintimes1, onbincount] = SliceHistogram(microtimes0, times0,
                                                                                                    limits0, microtimes1,
                                                                                                    times1, limits1, dtmicro,
                                                                                                    dtmacro, Ilims[binnr],
                                                                                                    Ilims[binnr + 1])
    
    
                #
                #
                # Calculate G2 for a time bin
                """
                if g2bool == True:
                    fig = plt.figure()
                    [g2tlist, g2ylist, _, _] = MakeG2(bintimes0, bintimes1, dtmicro, binmicrotimes0, binmicrotimes1, 8e-9, 220)
                    fig.tight_layout()
                    # plt.show()
    
                    if savebool == True:
                        fig.savefig(
                            lifetime_folder + namelist[measnr] + "_G2_" + str(binnr) + "_nrbins" + str(nrbins) + ".png",
                            dpi=300)
                    plt.close(fig)
                """
                #
                #
                #
    
                # calculate average lifetimes
                fig = plt.figure()
                plt.title(['binnr:', str(binnr)])
    
                if binnr == 0:  # use bin with lowest intensity to get background counts (idea: lowest intensity usually means fastest decay. chance is high to avoid build-up of photons)
                    [tauavelist[measnr, binnr], Alist[measnr, binnr], ybg] = GetLifetime(
                        np.append(binmicrotimes0, binmicrotimes1), dtmicro, dtmacro, dttau, -1, histbinmultiplier, 0)
    
                    if savebool == True:
                        try:
                            fig.savefig(
                                lifetime_folder + namelist[measnr] + '_lifetime' + str(binnr) + "_nrbins" + str(nrbins) + '.pdf',
                                bbox_inches='tight', format='pdf')
                            fig.savefig(
                                lifetime_folder + namelist[measnr] + '_lifetime' + str(binnr) + "_nrbins" + str(
                                    nrbins) + '.png',
                                bbox_inches='tight', dpi = 300)
                        except OSError as e:
                            print(e)
                            warnings.warn("Create the output/pdf/ folder")
                        plt.close(fig)
    
    
                    bgcts = ybg * int(dtmacro / (
                                dtmicro * histbinmultiplier))  # ybg is the background counts per histogram bin of lifetime binning. multiply with number of hist bins (int(dtmacro/(dtmicro*histbinmultiplier))) to get total nr of bgcounts
                    bgcps = bgcts / (onbincount * binwidth)  # divide by nr of bins of this slice times the binwidth to get cps
                    ybg = ybg / onbincount  # normalize to nr of time bins within slice
                    print('bgcts:', bgcts)
                    print('bgcps:', bgcps)
                    print('onbincount:', onbincount)
                    state_files.at[no_cols_df, "bgcts"] = bgcts
                    state_files.at[no_cols_df, "bgcps"] = bgcps
                    state_files.at[no_cols_df, "onbincount"] = onbincount
    
                    # convert data for fitting with max likelihood
                    bgcpsfit = bgcps
                    bgcpbfit = bgcps * binwidth / int(
                        dtmacro / (dtmicro * histbinmultiplier))  # counts per fitting bin (see below!)
                    print('bgcpbfit', bgcpbfit)
                else:  # for all other bins just extract lifetimes
                    [tauavelist[measnr, binnr], Alist[measnr, binnr], _] = GetLifetime(
                        np.append(binmicrotimes0, binmicrotimes1), dtmicro, dtmacro, dttau, -1, histbinmultiplier,
                        ybg * onbincount) #FIXME: check the calculation of the Background level again
                    fig.savefig(
                        lifetime_folder + namelist[measnr] + '_lifetime' + str(binnr) + "_nrbins" + str(nrbins) + '.pdf',
                        bbox_inches='tight', format='pdf')
                    fig.savefig(
                        lifetime_folder + namelist[measnr] + '_lifetime' + str(binnr) + "_nrbins" + str(nrbins) + '.png',
                        bbox_inches='tight', dpi = 300)
                    print('onbincount:', onbincount)
                    plt.close(fig)
    
            # get lifetimes for each time bin
            if singlebintaubool == True:
                MaxLikelihoodFunction_c = nb.jit(nopython=True)(MaxLikelihoodFunction)
                for tbinnr in tqdm(range(nrtbins - 1)):
                    microtimes = np.append(microtimes0[limits0[tbinnr]:limits0[tbinnr + 1]],
                                           microtimes1[limits1[tbinnr]:limits1[tbinnr + 1]])
                    [ylist, xlist] = np.histogram(microtimes, int(dtmacro / (dtmicro * histbinmultiplier)),
                                                  [0, int(dtmacro / dtmicro)])
                    tlist = (xlist[:-1] + 0.5 * (xlist[1] - xlist[0])) * dtmicro * 1e9  # convert x-axis to time in ns
                    idxstart, = np.where(tlist == next((x for x in tlist if x > tmax), 1))  # find index of dttau
                    istart = idxstart.item()  # convert index to normal scalar
                    idxend, = np.where(tlist == next( (x for x in tlist if x > dttau * 1e9), len(tlist)))  # find index of dttau
    
                    iend = idxend.item() + istart  # convert index to normal scalar
                    bgcpbfit = bgcps * binwidth / int(dtmacro / (dtmicro * histbinmultiplier))  # counts per fitting bin
                    try:
                        [tau, A] = MaxLikelihoodFit_c(tlist, ylist, istart, iend, bgcpbfit, MaxLikelihoodFunction_c, False)
                    except (ZeroDivisionError) as e:
                        print(e)
    
                    taubinlist[tbinnr] = tau
                    intbinlist[tbinnr] = len(microtimes)
    
                    # TODO: The Plotting does not what I would expect it to do
                    """ 
                    if np.mod(tbinnr,1000) == 0: #for every 1000th bin output a plot of the fit (just to visualize the process)
                        print('maxlikelihood tbinnr:',tbinnr)
                        [tau,A] = MaxLikelihoodFit_c(tlist,ylist,istart,iend,bgcpbfit, True) # do again just for plotting
                    """
    
                # write out lifetimes and intensities
                if savebool == True:
                    try:
                        np.save(folder + 'output/' + namelist[measnr] + '_FLID_taubin', taubinlist)
                        np.save(folder + 'output/' + namelist[measnr] + '_FLID_intbin', intbinlist)
                    except FileNotFoundError as e:
                        print(e)
                        warnings.warn("Manually create the 'output/' folder")
            #% PLOT DATA
            print("Starting Plotting...")
    
            # input
            savebool = True  # safe figures to folder /output/pdf/
            intmin = 0
            intmax = 1.1 * np.max(inttrace)
            intbins = 50
            taumin = 0
            taumax = 60 #nanoseconds
            taubins = 61
            nrmeas = len(namelist)
    
            # initialize FLID array
            #FLID = np.ndarray(shape=(intbins, taubins), dtype=int)  # init list to save intensity per bin
    
            #fig2 = plt.figure(figsize=(8, 4))
            #plt.xlabel('time (ns)')
            #plt.ylabel('normalized count rate')
    
            try:
                decaytlist = np.load(folder + 'output/' + namelist[measnr] + '_decay_tlist.npy')
                decayylist = np.load(folder + 'output/' + namelist[measnr] + '_decay_ylist.npy')
                #plt.figure()
                #plt.semilogy(decaytlist, decayylist / np.max(decayylist))
                #plt.title("Decay")
                #plt.show()
            except OSError as e:
                pass
                #print(e)
                #warnings.warn("The necessary numpy arrays have not been created. Set the 'singlebintaubool' variable to True")
    
            try:
                intbinlist = np.load(folder + 'output/' + namelist[measnr] + '_FLID_intbin.npy')
                taubinlist = np.load(folder + 'output/' + namelist[measnr] + '_FLID_taubin.npy')
            except OSError as e:
                pass
                #print(e)
                #warnings.warn("The necessary numpy arrays have not been created. Set the 'singlebintaubool' variable to True")
    
            timelist = binwidth * np.arange(len(intbinlist))
            """
            fig = plt.figure(figsize=(18, 10))
            gs = gridspec.GridSpec(2, 3, width_ratios=[1.5, 3, 1.5])
    
            # FLID MAP
            plt.subplot(gs[0])
            plt.title('', y=0.85)
            plt.ylabel('emission intensity (cts / %i ms)' % (binwidth * 1e3))
            plt.xlabel('PL lifetime (ns)')
            [H, intedges, tauedges, img] = plt.hist2d(taubinlist, intbinlist, [taubins, intbins],
                                                      [[taumin, taumax], [intmin, intmax]], cmin=0, cmap='RdPu',
                                                      norm=mpl.colors.LogNorm())
            plt.colorbar(img)
    
            # Intensity trace
            plt.subplot(gs[1])
            plt.plot(timelist, intbinlist, '-', linewidth=1.5, color='k')
            plt.xlabel('time (s)')
            plt.ylabel('counts / %i ms' % (binwidth * 1e3))
            plt.xlim([0, Texp])
            plt.ylim([intmin, intmax])
    
            # Intensity histogram
            inthistogram = np.histogram(intbinlist, intbins, [0, int(np.max(intbinlist))])
            plt.subplot(gs[2])
            xdata = 0.5 * (inthistogram[1][:-1] + inthistogram[1][1:])
            plt.plot(inthistogram[0], xdata, color='k', linewidth=1.5)
            plt.xlabel('occurrence')
            plt.ylabel('counts / %i ms' % (binwidth * 1e3))
            plt.ylim([intmin, intmax])
    
            # poisson fits
            #    popt = singlepoissonfit(xdata[int(intbins/2):-1],inthistogram[0][int(intbins/2):-1],[np.sqrt(0.9*np.max(intbinlist)),0.9*np.max(intbinlist)]);
            #    plt.plot(singlepoisson(xdata,*popt),xdata)
            #    popt = singlepoissonfit(xdata[0:int(intbins/2)],inthistogram[0][0:int(intbins/2)],[np.sqrt(0.12*np.max(intbinlist)),0.12*np.max(intbinlist)]);
            #    plt.plot(singlepoisson(xdata,*popt),xdata)
    
            # g2 correlation
            if g2bool == True:
                g2tlist = np.load(folder + 'output/' + namelist[measnr] + '_G2_tlist.npy')
                g2ylist = np.load(folder + 'output/' + namelist[measnr] + '_G2_ylist.npy')
                plt.subplot(gs[3])
                plt.plot(g2tlist, g2ylist / np.max(g2ylist), color = "k", linewidth=1.5)
                plt.xlabel('delay (ns)')
                plt.ylabel('occurence (a.u.)')
                plt.ylim([0, 1])
            else:
                plt.subplot(gs[3])
                plt.xlabel('time (ns)')
                plt.ylabel('normalized count rate')
                for measnr in range(nrmeas):
                    decaytlist = np.load(folder + 'output/' + namelist[measnr] + '_decay_tlist.npy')
                    decayylist = np.load(folder + 'output/' + namelist[measnr] + '_decay_ylist.npy')
                plt.semilogy(decaytlist, decayylist / np.max(decayylist), 'k')
    
            # Lifetime trace
            plt.subplot(gs[4])
            plt.plot(timelist, taubinlist, '-', linewidth=1.5, color='k')
            plt.xlabel('time (s)')
            plt.ylabel('lifetime (ns)')
            plt.xlim([0, Texp])
            plt.ylim([taumin, taumax])
    
            # lifetime histogram
            tauhistogram = np.histogram(taubinlist, taubins, [taumin, taumax])
            xdata = 0.5 * (tauhistogram[1][:-1] + tauhistogram[1][1:])
            plt.subplot(gs[5])
            plt.plot(tauhistogram[0], 0.5 * (tauhistogram[1][:-1] + tauhistogram[1][1:]), color='k', linewidth=1.5)
            plt.xlabel('occurrence')
            plt.ylabel('lifetime (ns)')
            plt.ylim([taumin, taumax])
    
            #poisson fit #FIXME: This does not work (it just produces a vertical line)
            #popt = singlepoissonfit(xdata,tauhistogram[0],[np.sqrt(0.95*np.max(taubinlist)),0.95*np.max(taubinlist)]);
            #plt.plot(singlepoisson(xdata,*popt),xdata)
    
            plt.suptitle(namelist[measnr])
            fig.tight_layout()
            if savebool == True:
                try:
                    print("saving...")
                    fig.savefig(folder + 'output/pdf/' + namelist[measnr] + '_2by3.pdf', bbox_inches='tight', format='pdf')
                except OSError as e:
                    print(e)
                    warnings.warn("Create the output/pdf/ folder")
    
            #plt.show()
            """
            print("")
    
        state_files.to_csv(g2_folder + "state_files.csv", index=False)
    print("All plots done")