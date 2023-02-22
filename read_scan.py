 # -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 13:52:43 2019
Reads data from Scanning_imaging and plots an image in the current form. Uncommenting in the for loop makes it also plot
lifetime images. The find particles section finds bright spots where the threshold should be significantly above the
noise floor and saves the index of their positions.

@author: rober
adapted for BSc Thesis by Julian Santen (from Sept. 2022)
"""
from functionDefinitions.functions_read_scan import *

if __name__ == "__main__":

    ########################################################
    # Define Name of Input Folder Here
    ########################################################
    folder = "C:/Users/omel-acton/Desktop/Gabriel_Local_EMCCD/22-11-18/"
    Mname= "AP-3-108_C1x10-6_SingleDot_1x10+6Hz_405exct_220nW_x_49_y78_z124_range40_pts200"

    MAX_SIGMA = 10 # ensure that a single blob is not identified as two individual blobs
    MIN_SIGMA = 1.5
    THRESHOLD = 200  # the higher, the fewer particles found
    savebool = True
    globalPlotBool = False  # if false *most* plots are not generated

    ########################################################
    # LEAVE UNCHANGED
    ########################################################

    name = Mname+'_tttrmode'
    HHsettings = load_obj(Mname+'_HHsettings', folder)
    MCLdata = load_obj(Mname+'_MCLdata', folder)
    xbinning = 1
    ybinning = 1

    # Check if Necessary Folders are created (namely: /output/pdf)
    output_folder = os.path.join(folder, 'output')
    pdf_folder = os.path.join(output_folder, 'pdf')
    png_folder = os.path.join(output_folder, 'png')

    if not (os.path.exists(output_folder)):  # if folder not found
        os.mkdir(output_folder)
        print("output folder created: ", output_folder)
    else:
        print("output folder found: ", output_folder)

    if not (os.path.exists(pdf_folder)):  # if folder not found
        os.mkdir(pdf_folder)
        print("pdf folder created: ", pdf_folder)
    else:
        print("pdf folder found: ", pdf_folder)

    if not (os.path.exists(png_folder)):  # if folder not found
        os.mkdir(png_folder)
        print("pdf folder created: ", png_folder)
    else:
        print("pdf folder found: ", png_folder)

    #
    # read data
    #
    [dtmicro, dtmacro, microtimes0, macrotimes0, microtimes1, macrotimes1, nrphotons0,nrphotons1,overflows,microtimesP,macrotimesP,nrP,microtimesL,macrotimesL,nrL] = importT3marked(folder + name + '.out', HHsettings) # import data from .ptu file
    Texp = round(overflows*dtmacro*1024,0) # macrotime values can store up to 2^10 = 1024 before overflow. overflows*1024*dtmacro gives experiment time [s]

    #
    # Plot and Process
    #

    # time gating (select only photons that arrived after tthreshold)
    ylist,tlist=GetHistogram(microtimes0,dtmicro,dtmacro,10)

    if globalPlotBool == True:
        fig = plt.figure()
        plt.plot(tlist*1e9,ylist,'.g')
        plt.xlabel('time (ns)')
        #plt.show()
        if savebool == True:
            fig.savefig(folder + "output/pdf/" + name + "_GetHistogram_all.pdf")
            fig.savefig(folder + "output/png/" + name + "_GetHistogram_all.png", dpi = 300, transparent = False)
        plt.close("all")


    tthreshold = 0*20000e-9 #threshold time in [s] #TODO: What is the meaning of this time?
    [__,__,microtimes0,macrotimes0] = TimeGate(microtimes0,macrotimes0,dtmicro,tthreshold)

    # get number of first photon per trigger pulse
    limits0P = FindFirstPhotonAfterTrigger(macrotimes0,macrotimesP)
    limits0L= FindFirstPhotonAfterTrigger(macrotimes0,macrotimesL)

    # Plot decay histogram for all photons
    if globalPlotBool == True:
        plt.figure(1)
        fitOut = GetLifetime(microtimes0,dtmicro,dtmacro,20e-9,plotbool=True,method='WLS',ybg=-1) # WLS = Weighted least squares fit


    # Numba booster (c compilation of maximum-likelihood function)
    if globalPlotBool == True:
        MaxLikelihoodFunction_c = nb.jit(nopython=True)(MaxLikelihoodFunction) # implement c function
        fitOut = GetLifetime(microtimes0,dtmicro,dtmacro,20e-9,plotbool=True,method='ML_c',ybg=-1) # compile c function for the first time


    # Get intensities per pixel
    ylen = int(nrL)
    xlen = int(nrP/nrL)
    #xlen = MCLdata["Points"]
    Imat = np.zeros((xlen,ylen))
    dwelltimes= np.zeros((xlen,ylen))
    for l in range(ylen):
        for p in range(xlen):
            if p>0: # this is a special case because there is no pixel trigger at start of each line, but instead a line trigger
                Imat[p,l] = limits0P[p+xlen*l]-limits0P[p+xlen*l-1]
                dwelltimes[p,l] = macrotimesP[p+xlen*l]-macrotimesP[p+xlen*l-1]
            else:
                Imat[p,l] = limits0P[p+xlen*l]-limits0L[l]
                dwelltimes[p,l] = macrotimesP[p+xlen*l]-macrotimesL[l]

    # Get lifetimes and intensities (binned)
    taumat = np.zeros((int(xlen/xbinning),int(ylen/ybinning)))
    Imat_bin = np.zeros((int(xlen/xbinning),int(ylen/ybinning)))
    Iclustmat_bin = np.zeros((int(xlen/xbinning),int(ylen/ybinning)))
    xcnt = 0;
    ycnt = 0;

    start_time = time.time()
    for l in range(0,ylen,ybinning):
        print('line',l,'out of',ylen-1) #start counting at 0
        for p in range(0,xlen,xbinning):
            if p>0: # this is a special case because there is no pixel trigger at start of each line, but instead a line trigger
                microtimesp = microtimes0[limits0P[p+xlen*l-1]:limits0P[p+(xbinning-1)+xlen*l]]
                Ip = limits0P[p+xlen*l]-limits0P[p+xlen*l-1]
                dwelltime=macrotimesP[p]-macrotimesP[p-1]
                for dl in range(1,ybinning-1):
                    np.append(microtimesp,microtimes0[limits0P[p+xlen*(l+dl)-1]:limits0P[p+(xbinning-1)+xlen*(l+dl)]])
                    Ip += limits0P[p+(xbinning-1)+xlen*(l+dl)] - limits0P[p+xlen*(l+dl)-1]
            else:
                microtimesp = microtimes0[limits0L[l]:limits0P[p+(xbinning-1)+xlen*l]]
                Ip = limits0P[p+(xbinning-1)+xlen*l]-limits0L[l]
                for dl in range(1,ybinning-1):
                    np.append(microtimesp,microtimes0[limits0L[l+dl]:limits0P[p+(xbinning-1)+xlen*(l+dl)]])
                    Ip += limits0P[p+(xbinning-1)+xlen*(l+dl)] - limits0L[l+dl]

            # TODO: What is the purpose of this figure? It slows down the code drastically
            """
            plt.figure()
            fitOut = GetLifetime(microtimesp,dtmicro,dtmacro,20e-9,plotbool=False,method='ML_c') # ML = Maximum likelihood fit
            taumat[xcnt,ycnt] = fitOut[0]
            """

            Imat_bin[xcnt,ycnt] = Ip
            Iclustmat_bin[xcnt,ycnt] = Ip>500
            xcnt += 1
        xcnt = 0 # reset column counter when reaching new line
        ycnt+=1
    print("--- %s s for execution of lifetime fitting ---" % (time.time() - start_time))

    # %% Timegated PL images not yet operational
    #Imat_bing = np.zeros((int(xlen/xbinning),int(ylen/ybinning)))
    #xcnt = 0;
    #ycnt = 0;
    #measstart=10*1e-9/dtmicro #start is at 10 ns
    #timemin=5*1e-9/dtmicro
    #timemax=10e-9/dtmicro
    #myphotons = np.array([i for i in range(len(microtimes0)) if (measstart+timemin)<microtimes0[i]<(measstart+timemax)])
    #for l in range(0,ylen,ybinning):
    #    print('line',l,'out of',ylen-1) #start counting at 0
    #    for p in range(0,xlen,xbinning):
    #        if p>0: # this is a special case because there is no pixel trigger at start of each line, but instead a line trigger
    #            microtimesp = myphotons[limits0P[p+xlen*l-1]:limits0P[p+(xbinning-1)+xlen*l]]
    #            Ipg = limits0P[p+xlen*l]-limits0P[p+xlen*l-1]
    #            for dl in range(1,ybinning-1):
    #                np.append(microtimesp,myphotons[limits0P[p+xlen*(l+dl)-1]:limits0P[p+(xbinning-1)+xlen*(l+dl)]])
    #                Ipg += limits0P[p+(xbinning-1)+xlen*(l+dl)] - limits0P[p+xlen*(l+dl)-1]
    #        else:
    #            microtimesp = myphotons[limits0L[l]:limits0P[p+(xbinning-1)+xlen*l]]
    #            Ipg = limits0P[p+(xbinning-1)+xlen*l]-limits0L[l]
    #            for dl in range(1,ybinning-1):
    #                np.append(microtimesp,myphotons[limits0L[l+dl]:limits0P[p+(xbinning-1)+xlen*(l+dl)]])
    #                Ipg += limits0P[p+(xbinning-1)+xlen*(l+dl)] - limits0L[l+dl]
    #        #plt.figure()
    #        Imat_bing[xcnt,ycnt] = Ipg
    #        xcnt += 1
    #    xcnt = 0 # reset column counter when reaching new line
    #    ycnt+=1
    #xvec=np.mean(MCLdata["xposmatrix"],1)
    #yvec=np.mean(MCLdata["yposmatrix"],0)
    ##plt.pcolor(xvec,yvec,np.transpose(Imat/(dwelltimes*dtmacro)))
    ##plt.pcolor(yvec,xvec,MCLdata["xpos"])
    #plt.imshow(Imat_bing)
    #plt.xlabel('y (um)')
    #plt.ylabel('x (um)')
    #plt.colorbar()


    #%% Plotting Definitions
    scancenter = [0,0]
    scancenter = MCLdata["ScanCenter"]
    scanrange = 20     #in micrometer %TODO: Is that a fixed value?
    #scanrange = MCLdata["ScanRange"]
    #points = 160
    points=MCLdata["Points"]
    xtargetvals=list(scancenter[0]+((j+0.5-points/2)/points)*scanrange for j in range(points)) #TODO: These are not used
    ytargetvals=list(scancenter[1]+((j+0.5-points/2)/points)*scanrange for j in range(points))
    xvec = np.mean(MCLdata["xposmatrix"],1)
    yvec = np.mean(MCLdata["yposmatrix"],0)


    #%% Plot intensity and lifetime maps
    if globalPlotBool == True:
        fig = plt.figure(figsize = (7,5))

        plt.subplot(221)
        plt.imshow(Imat_bin,extent=[min(xvec),max(xvec),min(yvec),max(yvec)])
        plt.colorbar()
        plt.title('PL intensity (counts per bin)')
        plt.xlabel('x (um)')
        plt.ylabel('y (um)')

        # Plot FLID map
        plt.subplot(222)
        plt.hist2d(Imat_bin.reshape(Imat_bin.size,1)[:,0],taumat.reshape(taumat.size,1)[:,0],(50,50))
        plt.title('FLID map')
        plt.xlabel('Counts per bin')
        plt.ylabel('Lifetime (ns)')
        plt.tight_layout()
        plt.colorbar()
        #plt.gca().set_aspect('equal', 'datalim')

        plt.subplot(223)
        plt.imshow(taumat,vmin=4,vmax=7,extent=[min(xvec),max(xvec),min(yvec),max(yvec)])
        plt.colorbar()
        plt.title('Lifetime (ns)')
        plt.xlabel('x (um)')
        plt.ylabel('y (um)')

        fig.tight_layout()
        plt.show()

        if savebool == True:
            #fig.tight_layout()
            fig.savefig(folder + "output/pdf/" + name + "_PL_intensity.pdf")
            fig.savefig(folder + "output/png/" + name + "_PL_intensity.png", dpi = 300, transparent = False)


    #%% plot triggers
    if globalPlotBool == True:
        fig = plt.figure(2)
        plt.plot(macrotimesP,np.ones(macrotimesP.shape),'r.')
        plt.plot(macrotimesL,np.ones(macrotimesL.shape),'b.')

        fig.tight_layout()
        plt.show()

        if savebool == True:
            fig.savefig(folder + "output/pdf/" + name + "_Triggers.pdf")
            fig.savefig(folder + "output/png/" + name + "_Triggers.png", dpi=300, transparent=False)

    testarray = np.array([macrotimesP[i+1]-macrotimesP[i] for i in range(len(macrotimesP)-1)])

    #%% Find particles
    #
    #


    #fig, ax = plt.subplots()
    #ax.imshow(Imat/(dwelltimes*dtmacro)) #FIXME: Why is the axis scaling here different to the pictures above?

    #a = np.hstack(Imat/(dwelltimes*dtmacro))
    #plt.hist(a, bins='auto')  # arguments are passed to np.histogram
    #plt.show()

    testimage = Imat/(dwelltimes*dtmacro)

    # https://scikit-image.org/docs/stable/api/skimage.feature.html#skimage.feature.blob_log
    # The blobs are returned Brightes -> Darkest
    # A 2d array with each row representing 2 coordinate values plus the sigma used.
    blobs_log = blob_log(testimage, min_sigma=MIN_SIGMA, max_sigma=MAX_SIGMA, threshold=THRESHOLD)  # log := laplacian of gaussians
    blobs_log[:, 2] = blobs_log[:, 2] * sqrt(2)  # calculate the approximate radius of the blob

    #%%
    print(str(blobs_log.shape[0])," Particles found")
    MCLdata["ParticleLocationsIndices"] = blobs_log
    save_obj(MCLdata, Mname+"_MCLdata", folder )

    # %% Plot PL Intensity Single Plot
    if True:
        fig, ax = plt.subplots()
        plt.imshow(Imat_bin, extent=[min(xvec), max(xvec), min(yvec), max(yvec)])
        plt.rcParams.update({'axes.titlesize': 'medium'})
        #plt.suptitle("PL intensity (counts per bin)")
        plt.title('Blobs: {} - Threshold: {} - Sigma $\in$ [{},{}]'.format(str(blobs_log.shape[0]), THRESHOLD, MIN_SIGMA, MAX_SIGMA))
        plt.xlabel('x [$\mu$m]')
        plt.ylabel('y [$\mu$m]')
        plt.gca().set_aspect('equal', 'datalim')
        plt.colorbar()

        # Add the circles for the detected particles (this is not validated, JPS)
        for j in range(MCLdata["ParticleLocationsIndices"].shape[0]):
            y, x, r = MCLdata["ParticleLocationsIndices"][j]
            y = int(y)
            x = int(x)

            # Based on logic. This works well, but should not... (code taken from above)
            xpos = x / (xvec.shape[0] - 1) * (xvec.max() - xvec.min()) + xvec.min()
            ypos = - y / (yvec.shape[0] - 1) * (yvec.max() - yvec.min()) + yvec.max()
            c = plt.Circle((xpos, ypos), r / 10, color='w', linewidth=1, fill=False,alpha=0.7)
            ax.add_patch(c)
            #plt.text(xpos, ypos, str(int(xpos)) + ";" + str(int(ypos)), color='white', alpha=0.5)  # add position

        fig.tight_layout()

        if savebool == True:
            fig.savefig(folder + "output/pdf/" + name + "_" + str(MAX_SIGMA) + "_" + str(THRESHOLD) + "_PL_Intensity_1by1.pdf")
            fig.savefig(folder + "output/png/" + name + "_" + str(MAX_SIGMA) + "_" + str(THRESHOLD) + "_PL_Intensity_1by1.png", dpi=300, transparent=False)

        plt.close("all")

    #%% Check some code from the measure_particles.py file
    """
    MCLdata = load_obj(Mname + '_MCLdata', folder)

    for j in range(MCLdata["ParticleLocationsIndices"].shape[0]):
        index0, index1, r = MCLdata["ParticleLocationsIndices"][j]
        index0 = int(index0)
        index1 = int(index1)
        print(np.round(MCLdata["xposmatrix"][index0, index1],2), np.round(MCLdata["yposmatrix"][index0, index1],2))
        """

    #%% Plot cross sections # FIXME: This does not work and I dont understand what's the intention
    #fig5 = plt.figure(5)
    #from scipy.optimize import curve_fit
    #from scipy import asarray as ar,exp
    #def gauss_function(x, a, x0, sigma):
    #    return a*np.exp(-(x-x0)**2/(2*sigma**2))
    #x = xvec
    #y = Imat[56,:]/(dwelltimes[56,:]*dtmacro)
    #plt.plot(x,y, "-*", color = "red")
    #n = len(x)                          #the number of data
    #mean = sum(x*y)/n                   #note this correction
    #sigma = sum(y*(x-mean)**2)/n        #note this correction
    #popt, pcov = curve_fit(gauss_function, x, y, p0 = [1000, 156, 0.2])
    #plt.plot(x, gauss_function(x, *popt), label='fit', color = 'blue')
    #plt.show()

    #%% Save arrays to hard drive
    if savebool == True:
        np.save(folder +'output/' + name +'_taumat',taumat)
        np.save(folder +'output/' + name +'_Imat',Imat)
        np.save(folder +'output/' + name + '_Imat_bin',Imat_bin)

    print("All Done")