"""
created by Julian Santen on 11.11.2022 for BSc

Purpose: Analyse spectra obtained in andor of potential single dots and fit  a gaussian or lorentzian to their emission.
The spectrum should be exported in andor to individual .asc files. It is assumed that each file contains only one
relevant quantum dot! (This is a key assumption! If this does not hold, a peak finding function would need to be
implemented)
"""
import os
import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt

matplotlib.use("Qt5Agg")
import sys
import warnings
import seaborn as sns
import time
import re

sns.set()
sns.color_palette("bright")
sns.set_color_codes(palette='bright')

#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#rc('text', usetex=True)

import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
from lmfit import Model
from lmfit.models import GaussianModel, ConstantModel, LorentzianModel, LinearModel, QuadraticModel, PolynomialModel

########################################################
# Function Definitions
########################################################


"""
Create a new dataframe (df). The df_fit dataframe contains only a narrow wavelength region that will be used for a 
gaussian / lorentzian fit. This region is constrained by [min_lambda_fit; max_lambda_fit].
The df_plot frame contains more data, to put the fit data perspective.
"""


def create_df(filepath: str, min_lambda_fit: float, max_lambda_fit: float):
    """
    Create a dataframe based on a provided csv file (that contains the andor data). Also integrate / sum the andor
    spectrum in the vertical direction. This assumes that only a single blob is scanned in each andor file
    :param filepath: filepath to specific csv file
    :param min_lambda_fit: minimum wave length that will be used for fit
    :param max_lambda_fit: maximum wave length that will be used for fit
    :return: return the dataframe with integrated values and calculated emission energy
    """
    df = pd.read_csv(filepath, header=None)
    df = df.drop(df.columns[-1:], axis=1)  # drop the last column because it contains NaNs

    # create a df that only contains a region that is relevant due to optical filters
    df_fit = df[(df.iloc[:, 0] > min_lambda_fit) & (df.iloc[:, 0] < max_lambda_fit)]

    # sum / integrate over the y direction
    df_fit["integral"] = (df_fit.iloc[:, 1:].values.sum(axis=1) / np.size(df_fit.iloc[:, 1:], axis=1))

    # add the emission energies (from: E = h * c / lambda)
    df_fit["EmissionEnergy"] = 4.135667696e-15 * 299792458 * 1e9 / df.iloc[:, 0]

    return df_fit





def create_folder(folderpath: str):
    """
    Create a folder given a path to a folder. This does not work recursively atm (i.e., you have to create every level of a
    nested folder step by step
    """
    if not (os.path.exists(folderpath)):  # if folder not found
        os.mkdir(folderpath)  # create path
        print("folder created: ", folderpath)
    return folderpath


def fit_peak_model(yval, xval, bool_gaussian, bool_lorentzian, bool_weight=True):
    """
    Fit a peak model (lorentzian or gaussian) to the input data
    :param yval: contains y data for the fit
    :param xval: contains x data for the fit
    :param bool_gaussian: if True, a Gaussian Curve will be fitted (while bool_lorentzian == False)
    :param bool_lorentzian:  if True, a Lorentzian Curve will be fitted (while bool_gaussian == False)
    :param bool_weight: if True, a weighted fit will be solved. The weight corresponds to the normalized y-value and
    therefore puts more emphasis on data that is further from the x-axis
    :return: returns the result of the fit. this file format is from the lmfit package
    """
    if bool_gaussian & bool_lorentzian:
        print("Choose either Gaussian or Lorentzian fit - not both!")
        model = GaussianModel()
    elif bool_gaussian:
        model = GaussianModel()
    elif bool_lorentzian:
        model = LorentzianModel()
    else:
        return

    # TODO: This is not the best handling of infinities
    if not np.any(np.isinf(yval)) and not np.any(np.isinf(xval)):
        pass
    elif np.any(np.isinf(yval)) or np.any(np.isinf(xval)):
        warnings.warn("Triggered Infinity. Returns zero fit for this plot (can be ignored)")
        # mock fit (this is non-sensical data just so I can provide a return value)
        xval = np.ones(np.size(yval))
        yval = np.ones(np.size(yval))
        pars = model.guess(yval, x=xval, amplitude=-0.5)
        result = model.fit(yval, pars, x=xval)
        return result

    pars = model.guess(yval, x=xval, amplitude=-0.5)

    if bool_weight:
        weight = np.abs(yval) / np.max(yval)
        result = model.fit(yval, pars, x=xval, weights=weight)
    else:
        result = model.fit(yval, pars, x=xval)

    return result



def fit_bg_model(yval, xval, degree_bg_fit=3):
    """
    fit a polynomial function on the background data
    """
    model = PolynomialModel(degree=degree_bg_fit)

    pars = model.guess(yval, x=xval, amplitude=-0.5)
    result = model.fit(yval, pars, x=xval)

    return result



def iterate_and_fit_all_files(folderpath: str, save_folderpath: str, background_fit, bool_gaussian, bool_lorentzian,
                              min_energy, max_energy, g2df, bool_weight,
                              specific_files={}, excluded_files=[], chisqr_threshold=5):
    """
    fit the dataframe with a gaussian: to do that, find the background value and subtract it, then fit the curve
    optional: normalize the plot
    """

    # create Dataframes
    df_fit_coeffs = pd.DataFrame()  # to store the fit coefficients
    array_series_x = []  # to store the actual x-values of the spectra
    array_series_y = []  # to store the actual y-values of the spectra
    array_series_state = []  # store the state of the spectra
    array_series_counter_spectrum = []  # store an increasing counter for a specific day
    array_series_included = []  # store the threshold based on chi-sqr
    array_series_chisqr = []  # store the actual chi-sqr value
    array_series_height = []
    array_series_a_normed = []
    array_series_mean = []
    array_series_sigma = []
    array_series_fwhm = []


    # filename_export
    filename_lst = ""
    plot_name = "Gaussian" if bool_gaussian else "Lorentzian"
    print("Plot Name: ", plot_name)

    #
    # Determine which files will be analyzed based on the provided arrays
    #

    # if exclusion array given but no inclusion array: subtract exclusion array from all files
    if (specific_files is iterate_and_fit_all_files.__defaults__[0]) and (
            excluded_files is not iterate_and_fit_all_files.__defaults__[1]):
        all_files = np.array([])
        for filename in sorted(os.listdir(folder)):
            if filename.endswith(".asc"):
                all_files = np.append(all_files, int(filename[-8:-4]))
        specific_files = np.array(list(set(all_files) - set(excluded_files))).astype(int)
        print("exclusion array provided: ", excluded_files)
        print("chose: ", specific_files)

    # if inclusion array given, but no exclusion array: only analyzed the included files
    elif (specific_files is not iterate_and_fit_all_files.__defaults__[0]) and (
            excluded_files is iterate_and_fit_all_files.__defaults__[1]):
        specific_files = np.array(specific_files)
        print("array provided: ", specific_files)

    # else choose all arrays
    else:
        specific_files = np.array([])
        for filename in sorted(os.listdir(folder)):
            if filename.endswith(".asc"):
                specific_files = np.append(specific_files, int(filename[-8:-4]))
        print("array not provided choose all files")
        # print("array not provided choose all files: ", specific_files)

    # iterate over all files that are .asc()
    for filename in sorted(os.listdir(folderpath)):
        if filename.endswith(".asc") and int(filename[-8:-4]) in specific_files:
            particle_number = int(filename[-8:-4])
            print(particle_number)
            filename_lst = filename  # save the filename to export it later

            # create a dataframe
            df = create_df(os.path.join(folder, filename), min_lambda_fit, max_lambda_fit)

            # subtract the background level
            df["integral"] -= background_fit.eval(x=df["EmissionEnergy"].values)

            # fit with provided model
            result_peak = fit_peak_model(yval=df["integral"].values, xval=df["EmissionEnergy"].values,
                                         bool_gaussian=bool_gaussian, bool_lorentzian=bool_lorentzian, bool_weight=bool_weight)

            # store results in dataframe
            df_fit_coeffs.at[particle_number, 'a'] = result_peak.params['amplitude'].value
            df_fit_coeffs.at[particle_number, 'mean'] = result_peak.params['center'].value
            df_fit_coeffs.at[particle_number, 'sigma'] = result_peak.params['sigma'].value
            df_fit_coeffs.at[particle_number, 'redchi'] = result_peak.redchi
            df_fit_coeffs.at[particle_number, 'chisqr'] = result_peak.chisqr
            df_fit_coeffs.at[particle_number, 'fwhm'] = result_peak.params['fwhm'].value
            df_fit_coeffs.at[particle_number, 'height'] = result_peak.params['height'].value


            # normalize by peak height
            df["integral_normed"] = df["integral"] * 1 / np.abs(result_peak.params['height'].value)

            # fit again now on normalized plots
            result_normed = fit_peak_model(yval=df["integral_normed"].values, xval=df["EmissionEnergy"].values,
                                           bool_gaussian=bool_gaussian, bool_lorentzian=bool_lorentzian)

            # store results in dataframe
            df_fit_coeffs.at[particle_number, 'a_normed'] = result_normed.params['amplitude'].value
            df_fit_coeffs.at[particle_number, 'mean_normed'] = result_normed.params['center'].value
            df_fit_coeffs.at[particle_number, 'sigma_normed'] = result_normed.params['sigma'].value
            df_fit_coeffs.at[particle_number, 'fwhm_normed'] = result_normed.params['fwhm'].value
            df_fit_coeffs.at[particle_number, 'height_normed'] = result_normed.params['height'].value
            df_fit_coeffs.at[particle_number, 'redchi_normed'] = result_normed.redchi
            df_fit_coeffs.at[particle_number, 'chisqr_normed'] = result_normed.chisqr

            try:
                # add the g2 state if the corresponding andor files exists
                df_fit_coeffs.at[particle_number, 'state_g2'] = g2df[g2df["particle_no"] == particle_number]["state"].values[0]
                df_fit_coeffs.at[particle_number, 'g2_tau_zero'] = g2df[g2df["particle_no"] == particle_number]["g2_tau_zero"].values[0]
            except (KeyError, IndexError) as e:
                #print("Key Error particle_number: ", particle_number, filename)
                df_fit_coeffs.at[particle_number, 'state_g2'] = "not_provided"

            # check if chisqr values satisfy threshold and other possible constraints
            if (df_fit_coeffs.at[particle_number, 'chisqr_normed'] < chisqr_threshold) and \
                    (0 <= df_fit_coeffs.at[particle_number, 'height_normed'] <= 1.2) and \
                    (df_fit_coeffs.at[particle_number, 'a_normed'] <= 2) and \
                    (min_energy <= df_fit_coeffs.at[particle_number, 'mean'] <= max_energy) and \
                    (0.027 <= df_fit_coeffs.at[particle_number, 'sigma'] <= 1):

                df_fit_coeffs.at[particle_number, 'included'] = True
            else:
                df_fit_coeffs.at[particle_number, 'included'] = False

            # store the series values in designated data frame (no matter if it is included or not)
            array_series_x.extend(df["EmissionEnergy"].values)
            # print(array_series_y)
            array_series_y.extend(df["integral"].values)
            array_series_state.extend(np.full((1, np.size(df["integral"].values)),
                                              df_fit_coeffs.at[particle_number, 'state_g2'])[0].tolist())
            array_series_counter_spectrum.extend(
                np.full((1, np.size(df["integral"].values)), particle_number)[0].tolist())
            array_series_included.extend(
                np.full((1, np.size(df["integral"].values)), df_fit_coeffs.at[particle_number, 'included'])[0].tolist())
            array_series_chisqr.extend(
                np.full((1, np.size(df["integral"].values)), df_fit_coeffs.at[particle_number, 'chisqr_normed'])[0].tolist())
            array_series_fwhm.extend(
                np.full((1, np.size(df["integral"].values)), df_fit_coeffs.at[particle_number, 'fwhm'])[0].tolist())

            # add the remaining coefficients
            array_series_height.extend(
                np.full((1, np.size(df["integral"].values)), df_fit_coeffs.at[particle_number, 'height_normed'])[0].tolist())
            array_series_a_normed.extend(
                np.full((1, np.size(df["integral"].values)), df_fit_coeffs.at[particle_number, 'a_normed'])[
                    0].tolist())
            array_series_mean.extend(
                np.full((1, np.size(df["integral"].values)), df_fit_coeffs.at[particle_number, 'mean'])[0].tolist())
            array_series_sigma.extend(
                np.full((1, np.size(df["integral"].values)), df_fit_coeffs.at[particle_number, 'sigma'])[0].tolist())


            # add the name
            df_fit_coeffs.at[particle_number, 'filename'] = filename_lst

            # plot the data TODO This is without the background offset!
            fig = plt.figure()
            sns.set_palette("bright")

            label_particle = "No: " + str(particle_number) + \
                             " chisqr: " + str(np.round(df_fit_coeffs.at[particle_number, 'chisqr_normed'], 2)) + \
                             " FWHM: " + str(np.round(df_fit_coeffs.at[particle_number, 'fwhm'],2)) + \
                             " included: " + str(df_fit_coeffs.at[particle_number, 'included'])

            plt.plot(df["EmissionEnergy"].values, df["integral"].values, 'o', ms=6, alpha=0.7, label=label_particle)
            #sns.scatterplot(x=df["EmissionEnergy"].values, y=df["integral"].values, label=label_particle, alpha=0.7)
            plt.title("Emission of Particle " + str(particle_number))
            plt.plot(df["EmissionEnergy"].values, result_peak.best_fit, 'k--', label="fit")
            # plt.plot(df["EmissionEnergy"].values, result_peak_unweighted.best_fit, 'k--', label = "unweighted")
            plt.xlabel("Emission Energy [eV]")
            plt.ylabel("counts [a.u.]")
            plt.legend(loc=1)
            fig.tight_layout()
            fig.savefig(save_folderpath + "individual/" + filename[:-9] + "_curve_fit_" + plot_name + str(particle_number) + ".pdf")
            #fig.savefig(save_folderpath + "individual/" + filename[:-9] + "_curve_fit_" + plot_name + str(
            #    particle_number) + ".png", dpi=300)
            plt.close("all")

    # print("array_series_state: ", array_series_state)

    df_series = pd.DataFrame({'x': array_series_x, 'y': array_series_y, 'state': array_series_state,
                              'counter_spectrum': array_series_counter_spectrum, 'included': array_series_included,
                              'chisqr_normed': array_series_chisqr, 'height_normed': array_series_height,
                              'a_normed': array_series_a_normed, 'mean': array_series_mean, 'sigma': array_series_sigma,
                              'fwhm': array_series_fwhm})

    df_fit_coeffs.sort_index(ascending=True, inplace=True)  # sort with ascending particle number
    return df_fit_coeffs, filename_lst[:-9] + plot_name, df_series


def bg_level_not_normalized(folderpath: str, save_folderpath: str, min_lambda_fit: float, max_lambda_fit: float,
                            degree_bg_fit=3, save_individual=True, bg_threshold=800):
    # create dataframe for background fit, if necessary (needs to be outside the for loop)
    df_bg = pd.DataFrame()
    number_of_bg = 0  # number of bg fits performed
    filename_global = ""
    bg_array = []

    # iterate over all files that are .asc()
    for filename in sorted(os.listdir(folderpath)):
        if filename.endswith(".asc"):
            particle_number = int(filename[-8:-4])  # TODO: andor counting differs by 1 from python
            print(particle_number)
            filename_global = filename

            # create a dataframe
            df = create_df(os.path.join(folder, filename), min_lambda_fit, max_lambda_fit)

            if np.max(df["integral"]) > bg_threshold:  # skip information that is non-bg
                continue

            if number_of_bg == 0:
                df_bg["EmissionEnergy"] = df["EmissionEnergy"]
                df_bg["bg"] = df["integral"]
                number_of_bg += 1
                bg_array.append(particle_number)
            else:
                df_bg["bg"] = df_bg["bg"] + df["integral"]
                number_of_bg += 1
                bg_array.append(particle_number)

            if save_individual:
                # plot the data (add the offset )
                fig = plt.figure()
                sns.set_palette("bright")
                plt.plot(df["EmissionEnergy"].values, df["integral"].values, 'o', ms=6, alpha=0.7,
                         label=str(particle_number))
                plt.title("Emission of Particle " + str(particle_number))
                plt.xlabel("Emission Energy [eV]")
                plt.ylabel("counts [a.u.]")
                fig.tight_layout()
                fig.savefig(save_folderpath + "individual/" + filename[:-9] + "_used_for_bg_" + str(particle_number) + ".pdf")
                #fig.savefig(
                #    save_folderpath + "individual/" + filename[:-9] + "_used_for_bg_" + str(particle_number) + ".png",
                #    dpi=300)
                plt.close("all")

    try:
        df_bg["bg"] = df_bg["bg"] * 1 / number_of_bg
    except KeyError as e:
        print(folderpath)
        warnings.warn("The Background Threshold 'bg_threshold' is too low! Set it higher.")
        sys.exit()
    result_bg = fit_bg_model(yval=df_bg["bg"].values, xval=df_bg["EmissionEnergy"].values,
                             degree_bg_fit=degree_bg_fit)

    # plot the bg data (add the offset )
    fig = plt.figure()
    sns.set_palette("bright")
    plt.plot(df_bg["EmissionEnergy"].values, df_bg["bg"].values, 'o', ms=6, alpha=0.7)
    plt.title("Background Emission for " + str(number_of_bg) + " particles")
    plt.plot(df_bg["EmissionEnergy"].values, result_bg.eval(x=df_bg["EmissionEnergy"].values), "k--")

    plt.xlabel("eV")
    plt.ylabel("counts")
    fig.tight_layout()
    fig.savefig(save_folderpath + filename_global[:-9] + "_bg_fit" + ".pdf")
    plt.close("all")
    return result_bg, bg_array


def create_boxplot(df: pd.DataFrame, savepath: str, filename: str, bool_fit_gaussian):
    # https://towardsdatascience.com/scattered-boxplots-graphing-experimental-results-with-matplotlib-seaborn-and-pandas-81f9fa8a1801

    # Set style options here #
    sns.set_style("whitegrid")  # "white","dark","darkgrid","ticks"
    medianprops = dict(linewidth=1.5, linestyle='-', color='k')

    # set Plot Name
    plot_name = "Gaussian" if bool_fit_gaussian else "Lorentzian"

    # Make Sure Only Positive Sigmas
    df = df[df["included"] == True]

    fig, axes = plt.subplots(1, 4, figsize=(9.5 * 4 / 3, 4.2))
    sns.set_palette("bright")
    #rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    #rc('text', usetex=True)

    fig.subplots_adjust(top=0.3)

    plt.subplot(141)
    label1 = str(np.round(np.mean(df["a"]), 2)) + r" $\pm$ " + str(np.round(np.std(df["a"]), 2))

    plt.boxplot(df[["a"]].values, labels=[label1], medianprops=medianprops, showfliers=False)
    xs = np.random.normal(1, 0.04, df[["a"]].values.shape[0])  # adds jitter to the data points - can be adjusted
    plt.scatter(xs, df["a"], alpha=0.4, color="c")
    plt.ylabel("amplitude [a.u.]")

    plt.subplot(142)
    label2 = str(np.round(np.mean(df["mean"]), 2)) + " $\pm$ " + str(np.round(np.std(df["mean"]), 2)) + " eV"

    plt.boxplot(df[["mean"]].values, labels=[label2], medianprops=medianprops, showfliers=False)
    xs = np.random.normal(1, 0.04,
                          df[["mean"]].values.shape[0])  # adds jitter to the data points - can be adjusted
    plt.scatter(xs, df["mean"], alpha=0.4, color="r")
    plt.ylabel("center [eV]")

    plt.subplot(143)
    label3 = str(np.round(np.mean(df["sigma"]), 2)) + " $\pm$ " + str(
        np.round(np.std(df["sigma"]), 2)) + " eV"

    plt.boxplot(df[["sigma"]].values, labels=[label3], medianprops=medianprops, showfliers=False)
    xs = np.random.normal(1, 0.04,
                          df[["sigma"]].values.shape[0])  # adds jitter to the data points - can be adjusted
    plt.scatter(xs, df["sigma"], alpha=0.4, color="g")
    plt.ylabel("sigma [eV]")

    plt.subplot(144)
    label3 = str(np.round(np.mean(df["fwhm"]), 2)) + " $\pm$ " + str(
        np.round(np.std(df["fwhm"]), 2)) + " eV"

    plt.boxplot(df[["fwhm"]].values, labels=[label3], medianprops=medianprops, showfliers=False)
    xs = np.random.normal(1, 0.04,
                          df[["fwhm"]].values.shape[0])  # adds jitter to the data points - can be adjusted
    plt.scatter(xs, df["fwhm"], alpha=0.4, color="m")
    plt.ylabel("fwhm [eV]")

    plt.suptitle(
        plot_name + " Fit " + str(np.shape(df)[0]) + " Particles - " + filename)
    fig.tight_layout()

    fig.savefig(savepath + filename + "_curve_fit_" + plot_name + ".pdf")
    #fig.savefig(savepath + filename + "_curve_fit_" + plot_name + ".png", dpi=300)

    #
    # Also create correlation plots
    #

    sns.set_style("darkgrid")  # "white","dark","darkgrid","ticks"

    fig, axes = plt.subplots(1, 4, figsize=(9.5 * 4 / 3, 4.2))
    #rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    #rc('text', usetex=True)
    sns.set_palette("bright")

    fig.subplots_adjust(top=0.3)

    plt.subplot(141)
    plt.plot(df[["mean"]].values, df[["a"]].values, "o", color="c", alpha=0.3)
    plt.xlabel("center [eV]")
    plt.ylabel("amplitude [a.u.]")

    plt.subplot(142)
    plt.plot(df[["mean"]].values, df[["sigma"]].values, "o", color="g", alpha=0.3)
    plt.xlabel("center [eV]")
    plt.ylabel("sigma [eV]")

    plt.subplot(143)
    plt.plot(df[["mean"]].values, df[["fwhm"]].values, "o", color="m", alpha=0.3)
    plt.xlabel("center [eV]")
    plt.ylabel("fwhm [eV]")

    plt.subplot(144)
    plt.plot(df[["fwhm"]].values, df[["a"]].values, "o", color="r", alpha=0.3)
    plt.xlabel("fwhm [eV]")
    plt.ylabel("amplitude [a.u.]")

    plt.suptitle(
        "Correlation " + plot_name + " " + str(np.shape(df)[0]) + " Particles - " + filename)
    fig.tight_layout()
    fig.savefig(savepath + filename + "_correlation_" + plot_name + ".pdf")
    #fig.savefig(savepath + filename + "_correlation_" + plot_name + ".png", dpi=300)


# Compare the singular, non-singular, discarded plots from the general import
# def compare_singular_nonsingular():


# Class that handles the analysis of the files. This includes storing the the import files but also the results of
# of the Gaussian / Lorentzian fit as well as the average of the spectra
class SourceFiles:
    def __init__(self):
        self.folders = []  # contains the folders that *will be analyzed*
        self.folders_only_csv = []  # contains the folders for which the previously exported CSV file will be loaded
        self.min_energies = []
        self.max_energies = []
        self.bg_thresholds = []
        self.min_lambdas = []
        self.max_lambdas = []
        self.df = pd.DataFrame()
        self.df_series = pd.DataFrame()
        self.df_series_x = pd.DataFrame()  # contains the grouped series values in every column
        self.df_series_y = pd.DataFrame()
        # based on "state_files.csv" create arrays that contain the filenumber of {singular, non_singular_discarded}
        self.df_g2_files = pd.DataFrame()
        self.number_of_imports = 0  # count the number of times a file has been added to the class
        self.imports_to_be_analyzed = []  # contains the count number for files that will be analyzed

    def __repr__(self):
        return "Folders: " + str(self.folders) + "\n" + "min_energies: " + str(
            self.min_energies) + "\n" + "max_energies: " + str(self.max_energies) + "\n" + "bg_thresholds: " + str(
            self.bg_thresholds)

    def add_entry(self, folder, min_energy, max_energy, bg_threshold):
        """
        :param folder: folder where the files are stored including forward slash at the end
        :param min_energy: mininimum energy that should be used for a fit
        :param max_energy: maximum energy that should be used for a fit
        :param bg_threshold: if the maximum intensity is below this threshold value, it is assumed that the plot only
        contains background and will be used for the background fit
        :return:
        """
        # increment counter for number of imports
        self.number_of_imports += 1

        parent_folder = os.path.dirname(os.path.dirname(folder))  # get parent directory of folder

        # Check if g2_files (state_files.csv) exists in the respective g2 folder
        temp_path = os.path.join(parent_folder, "analysis", "g2") + "/state_files.csv"

        if os.path.isfile(temp_path):
            print("g2 files from 4_analyze_tcspc_2d_scan.py script found for: ", folder)

            # create temp dataframe and append that to the existing df at the end
            df_temp = pd.read_csv(temp_path)

            # extract the particle number from string name: regular expression to match the digits after _Particle_ddd_
            # ATTENTION: I add one because andor starts counting at zero and python starts counting at 0
            df_temp['particle_no'] = df_temp['name'].apply(
                lambda f: int(re.search(r"Particle_(\d{1,3})_", f).group(1)) + 1)

            df_temp['import_no'] = self.number_of_imports

            # append the g2 files to the existing dataframe
            self.df_g2_files = pd.concat([self.df_g2_files, df_temp], ignore_index=True)

            #print(self.df_g2_files.columns)
            #print(self.df_g2_files[['import_no', 'state', 'particle_no', 'name']])

        else:
            warn_text = "g2 files from 4_analyze_tcspc_2d_scan.py script *not* found for: " + folder
            warnings.warn(warn_text)

        # Check if the file already has been analyzed before using this scrip.
        # If it has been analyzed load the csv else analyze it # TODO: explanation is wrong
        temp_path = os.path.join(parent_folder, "analysis", "fit", "export")
        if os.path.isdir(temp_path):
            print("Spectra analysis files found for ", folder)
            self.folders_only_csv.append(
                folder)  # add the folder name, s.t. the previously exported CSV file will be loaded

            # create local output folders
            create_folder(os.path.dirname(os.path.dirname(folder)) + "/analysis/")
            create_folder(os.path.dirname(os.path.dirname(folder)) + "/analysis/fit/")
            create_folder(temp_path)

            # add the files that will be analyzed to the folder
            self.folders.append(folder)
            self.min_energies.append(min_energy)
            self.max_energies.append(max_energy)
            self.bg_thresholds.append(bg_threshold)
            self.min_lambdas.append(1240 / max_energy)
            self.max_lambdas.append(1240 / min_energy)
            self.imports_to_be_analyzed.append(self.number_of_imports)

        else:  # folder has *not* been analyzed before --> add it to the list
            print("Spectra analysis files not found for ", folder)

            # create local output folders
            create_folder(os.path.dirname(os.path.dirname(folder)) + "/analysis/")
            create_folder(os.path.dirname(os.path.dirname(folder)) + "/analysis/fit/")
            create_folder(temp_path)

            # add the files that will be analyzed to the folder
            self.folders.append(folder)
            self.min_energies.append(min_energy)
            self.max_energies.append(max_energy)
            self.bg_thresholds.append(bg_threshold)
            self.min_lambdas.append(1240 / max_energy)
            self.max_lambdas.append(1240 / min_energy)
            self.imports_to_be_analyzed.append(self.number_of_imports)
            print("create_folder for analysis scripts")

    def concat_df(self, df_call):
        """
        :param df_call: dataframe that should be appended to self.df
        """
        self.df = pd.concat([self.df, df_call], ignore_index=True)

    def average_spectra(self, filename: str, savepath: str, df=pd.DataFrame()):

        # https://stackoverflow.com/questions/51258954/python-3-6-get-average-y-for-all-same-x-coordinates
        sns.set_style("darkgrid")

        # if only df is given, this df contains a series with x and y values that need to be grouped
        # the grouped result is then concatenated to the member variables
        # This is the LOCAL call
        if df is not self.average_spectra.__defaults__[0]:  # if df *is* given (aka df is *not* the standard value)
            print("df provided in method")

            print("df size before group: ", np.size(df))
            df_group = df.groupby("x", as_index=False)["y"].mean()  # group df here
            print("df size after group: ", np.size(df_group))

            # add the grouped values (a) by vertical concatenating to the df_series dataframe
            self.df_series = pd.concat([self.df_series, df_group], ignore_index=True)

            # add the grouped values (b) by horizontal addition as rows to the df_series_x and df_series_y dataframes
            self.df_series_x = pd.concat([self.df_series_x, df_group["x"]], axis=1, ignore_index=True)
            self.df_series_y = pd.concat([self.df_series_y, df_group["y"]], axis=1, ignore_index=True)

            # retain backwards compatibility -- create new plots based on states
            df_group_singular = df[df["state"] == "singular"].groupby("x", as_index=False)["y"].mean()
            df_group_non_singular = df[df["state"] == "non_singular"].groupby("x", as_index=False)["y"].mean()
            df_group_discarded = df[df["state"] == "discarded"].groupby("x", as_index=False)["y"].mean()

            fig = plt.figure()
            #rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
            #rc('text', usetex=True)
            sns.set_palette("bright")
            plt.plot(df_group["x"], df_group["y"], "ko", alpha=0.1, label="all")
            plt.plot(df_group_singular["x"], df_group_singular["y"], "o", alpha=0.5, label="singular")
            plt.plot(df_group_non_singular["x"], df_group_non_singular["y"], "o", alpha=0.5, label="non-singular")
            plt.plot(df_group_discarded["x"], df_group_discarded["y"], "o", alpha=0.5, label="discarded")
            plt.xlabel("Emission Energy [eV]")
            plt.ylabel("counts [a.u.]")
            plt.legend()
            fig.tight_layout()
            fig.savefig(savepath + filename + "_local_average_spectra" + ".pdf")
            plt.close("all")

        else:
            # if df is not provided, the member dataframes df_series_x and df_series_y contain the series information in
            # each column (i.e., the first column of df_x and df_y belong together)

            # ROUND the already grouped plots to create an average
            self.df_series["x"] = np.round(self.df_series["x"].copy(), 2)
            self.df_series = self.df_series.copy().groupby("x", as_index=False)["y"].mean()

            fig = plt.figure()
            sns.set_palette("bright")

            # Also plot the individual particles
            for i in range(len(self.df_series_x.columns)):
                x = self.df_series_x.iloc[:, i].values
                y = self.df_series_y.iloc[:, i].values
                plt.plot(x, y, "o", label=str(i), alpha=0.3)

            plt.plot(self.df_series["x"], self.df_series["y"], "ko-", alpha=1, label="averaged")

            plt.xlabel("Emission Energy [eV]")
            plt.ylabel("counts [a.u.]")
            plt.legend()
            fig.tight_layout()
            fig.savefig(savepath + filename + "_global_average_spectra" + ".pdf")
            #fig.savefig(savepath + filename + "_global_average_spectra" + ".png", dpi=300)
            plt.close("all")

    def export_series_to_csv(self, savepath: str, filename: str):
        self.df_series.to_csv(savepath + filename + "_series_global.csv")
        self.df_series_x.to_csv(savepath + filename + "_series_x_global.csv")
        self.df_series_y.to_csv(savepath + filename + "_series_y_global.csv")
        self.df.to_csv(savepath + filename + "_0fit_coeffs_global.csv")

    def export_all_csv(self, savepath: str, filename: str, g2df: pd.DataFrame()):
        print("export path", savepath + "export/" + filename + "_series_global.csv")
        self.df_series.to_csv(savepath + "export/" + filename + "_series_global.csv")
        self.df_series_x.to_csv(savepath + "export/" + filename + "_series_x_global.csv")
        self.df_series_y.to_csv(savepath + "export/" + filename + "_series_y_global.csv")
        self.df_g2_files.to_csv(savepath + "export/" + filename + "_g2_states_GLOBAL.csv")  # this export might be empty
        g2df.to_csv(savepath + "export/" + filename + "_g2_states.csv")  # this export might be empty

########################################################
# Change Paramters
########################################################
bool_fit_gaussian = False
bool_fit_lorentzian = True
bool_create_indiv_plots = True
gaussian_chisqr_THRESHOLD = 5
lorentzian_chisqr_THRESHOLD = 5

if __name__ == "__main__":
    ########################################################
    # Import Files Here
    # First folder-entry is used for export of global plots
    ########################################################
    Import = SourceFiles()

    Import.add_entry(
        folder="C:/DataAnalysis/AP_3_108/22-11-08_AP-3-108_MultipleDots_2Batches/andor_export/",
        min_energy=2.25, max_energy=2.44, bg_threshold=1200)
    Import.add_entry(
        folder="C:/DataAnalysis/AP_3_108/22-11-11-AP-3-108_MultipleSingleDot_Measurement/andor_export/",
        min_energy=2.23, max_energy=2.46, bg_threshold=1000)
    Import.add_entry(
        folder="C:/DataAnalysis/AP_3_108/22-11-18-AP-3-108_MultipleSingleDot_Measurement/andor_export/",
        min_energy=2.23, max_energy=2.44, bg_threshold=1000)
    Import.add_entry(
        folder="C:/DataAnalysis/AP_3_108/22-11-23-AP_3_108/andor_export_2/",
        min_energy=2.17, max_energy=2.44, bg_threshold=1000)
    
    # Import.add_entry(
    #     folder="C:/Users/gabri/Documents/GitHub/Scanning_Fitting_GN_persona_rep/example_data/andor_export/",
    #     min_energy=2.171, max_energy=2.44, bg_threshold=800)


    # Files to be analyzed
    print("Import.imports_to_be_analyzed): ", Import.imports_to_be_analyzed)

    ########################################################
    # Leave Unchanged
    ########################################################
    print(Import)
    start_time = time.time()

    for folder, min_lambda_fit, max_lambda_fit, bg_threshold, min_energy, max_energy, import_no in \
            zip(Import.folders, Import.min_lambdas, Import.max_lambdas, Import.bg_thresholds,
                Import.min_energies, Import.max_energies, Import.imports_to_be_analyzed):
        print(folder)
        print("Import Number: ", import_no)

        # create local output folders
        _ = create_folder(os.path.dirname(os.path.dirname(folder)) + "/analysis/")
        save_folder_path_local = create_folder(os.path.dirname(os.path.dirname(folder)) + "/analysis/fit/")
        _ = create_folder(os.path.dirname(os.path.dirname(folder)) + "/analysis/fit/individual/")

        # find all .asc files in a given folder
        result_bg, array_bg = bg_level_not_normalized(folder, save_folder_path_local, min_lambda_fit, max_lambda_fit,
                                                      degree_bg_fit=3, bg_threshold=bg_threshold)

        # Check if an analysis for the g2 function has been found. If not, pass an empty dataframe to the function
        if Import.df_g2_files.empty:
            print("No g2 files for: ", folder)
            g2df = pd.DataFrame()
        else:
            g2df = Import.df_g2_files[Import.df_g2_files['import_no'] == import_no]

        df_fit_coefficients, filename_local, df_local_series = \
            iterate_and_fit_all_files(folder, save_folder_path_local,
                                      result_bg,
                                      chisqr_threshold=gaussian_chisqr_THRESHOLD,
                                      bool_gaussian=bool_fit_gaussian,
                                      bool_lorentzian=bool_fit_lorentzian,
                                      min_energy=min_energy,
                                      max_energy=max_energy,
                                      excluded_files=array_bg,
                                      bool_weight=False,
                                      g2df=g2df)

        # export the individual fits
        df_fit_coefficients.to_csv(save_folder_path_local + "export/" + filename_local + "_fit_coeffs_local.csv")
        df_local_series.to_csv(save_folder_path_local + "export/" + filename_local + "_series_local.csv", index=False)

        # create local boxplot
        create_boxplot(df_fit_coefficients, savepath=save_folder_path_local, filename=filename_local,
                       bool_fit_gaussian=bool_fit_gaussian)

        # create local average of spectra
        Import.average_spectra(df=df_local_series, filename=filename_local, savepath=save_folder_path_local)

        # add the fit coefficients to the overall dataframe
        Import.concat_df(df_fit_coefficients)

        # export analyzed data for easier access in future
        Import.export_all_csv(filename=filename_local, savepath=save_folder_path_local,
                              g2df=g2df)


    print("all done")
    print("--- %.2f seconds ---" % (time.time() - start_time))