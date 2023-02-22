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

# IGNORE: Matplotlibs MatplotlibDeprecationWarning (caused by shared y axis in side by side boxplot)
import warnings
warnings.filterwarnings("ignore",category=matplotlib.MatplotlibDeprecationWarning)

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np

class csvFiles:
    df = pd.DataFrame()  # overall DataFrame that contains the spectra data with some added attributes
    df_fit_coeffs = pd.DataFrame()  # contains the fit paramters from the previous Gaussian/Lorentzian fit

    parent_folder = ""
    timestr = time.strftime("%Y%m%d-%H%M%S")  # create a timestap with the current execution time

    def __init__(self):
        self.save_folder = save_folder_path_global

    def __repr__(self):
        pass

    def add_entry(self, group_number: int, group_name: str, folder: str, state="", color="", threshold_chisqr=None,
                  plot_individual_spectra=False, linestyle="", selected_blobs=[], alpha=2, upper_limit_fwhm=0.2):
        self.parent_folder = os.path.dirname(os.path.dirname(folder))  # get parent directory of folder

        """
        If no state is provided all files will be imported 
        
        Re-Import the previously created CSVs

        Make sure to:
        * BUGFIX: import number is currently wrong (check the exported data frame).
        * Think about rewriting the entire thing
        * include the fit coeffients
        * include the spectra information WITH information about the g2 fit. This might require a different g2 export
        """

        df_temp = pd.DataFrame()

        # create temporary df that will be appended to df
        self.parent_folder = os.path.join(os.path.dirname(os.path.dirname(folder)), "analysis", "fit", "export/")

        # import spectra
        for filename in os.listdir(self.parent_folder):

            # boxplots (aka fit coefficients are imported)
            if filename.endswith("fit_coeffs_local.csv") and "Lorentzian" in filename:
                df_temp_fit = pd.read_csv(self.parent_folder + filename)

                # extract the particle numbers
                df_temp_fit['particle_no'] = df_temp_fit['filename'].apply(
                    lambda f: int(re.search(r"(\d{1,4}).asc", f).group(1)))

                if state is self.add_entry.__defaults__[0]:  # if no state is provided choose all
                    pass

                elif state is not self.add_entry.__defaults__[0] and state in ["singular", "non_singular", "discarded"]:
                    df_temp_fit = df_temp_fit[df_temp_fit["state_g2"] == state]  # only select the state, chosen in the function call

                elif state is not self.add_entry.__defaults__[0] and state in [True, False]:# included state is given (this is a bool)
                    df_temp_fit = df_temp_fit[df_temp_fit["included"] == state]

                else:
                    warnings.warn(
                        'State not recognized. Either choose no state at all or "singular", "non_singular", \
                        "discarded" or the boolean values True or False (not a string!)')

                if selected_blobs is not self.add_entry.__defaults__[5]:  # if selected_blobs is provided
                    df_temp_fit = df_temp_fit[df_temp_fit['particle_no'].isin(selected_blobs)]  # select specified blobs

                if threshold_chisqr is not self.add_entry.__defaults__[2]:  # if chi_squared is provided
                    df_temp_fit = df_temp_fit[(df_temp_fit["chisqr_normed"] < threshold_chisqr) &
                                      (df_temp_fit["height_normed"] <= 1.2) &
                                      (df_temp_fit["a_normed"] <= 2) &
                                      (df_temp_fit["mean"] <= 2.6) &  # ATTENTION: Hardvoded energy cut-off
                                      (df_temp_fit["mean"] >= 2.0) &
                                      (df_temp_fit["sigma"] >= 0.027) &
                                      (df_temp_fit["sigma"] <= 1) &
                                      (df_temp_fit["fwhm"] <= upper_limit_fwhm)]

                    if np.shape(df_temp_fit)[0] == 0:
                        print("chisqr too low")


                df_temp_fit["group_number"] = group_number
                df_temp_fit["group_name"] = group_name

            elif filename.endswith("_series_local.csv") and "Lorentzian" in filename:  # find file that ends with _series_local #TODO: Make this part of the function call
                df_temp = pd.read_csv(self.parent_folder + filename)

                if state is self.add_entry.__defaults__[0]:  # if no state is provided choose all
                    #print("no state")
                    pass

                elif state is not self.add_entry.__defaults__[0]:  # if state is provided check what kind of state

                    if state in ["singular", "non_singular", "discarded"]:  # check if a g2 state is given
                        print("g2 state")
                        df_temp = df_temp[df_temp["state"] == state]  # only select the state, chosen in the function call

                    elif state in [True, False]:  # included state is given (this is a bool)
                        print("boolean state")
                        df_temp = df_temp[df_temp["included"] == state]

                else:
                    warnings.warn(
                        'State not recognized. Either choose no state at all or "singular", "non_singular", \
                        "discarded" or the boolean values True or False (not a string!)')

                if selected_blobs is not self.add_entry.__defaults__[5]:  # if selected_blobs is provided
                    # extract the particle numbers
                    print("selected blobs")
                    df_temp = df_temp[df_temp['counter_spectrum'].isin(selected_blobs)]  # select specified blobs

                if threshold_chisqr is not self.add_entry.__defaults__[2]:  # if chi_squared is provided
                    df_temp = df_temp[(df_temp["chisqr_normed"] < threshold_chisqr) &
                                      (df_temp["height_normed"] <= 1.2) &
                                      (df_temp["a_normed"] <= 2) &
                                      (df_temp["mean"] <= 2.6) &  # ATTENTION: Hardcoded energy cut-off
                                      (df_temp["mean"] >= 2.0) &
                                      (df_temp["sigma"] >= 0.027) &
                                      (df_temp["sigma"] <= 1) &
                                      (df_temp["fwhm"] <= upper_limit_fwhm)]

                df_temp["group_number"] = group_number
                df_temp["group_name"] = group_name

        # after iterating over all files
        if color is not self.add_entry.__defaults__[1]:  # if color is provided add this color
            df_temp["color"] = color  # add the color for the spectrum plot, if color is wanted

        if plot_individual_spectra is not self.add_entry.__defaults__[3]:  # if plot_individual_spectra provided: add it
            df_temp["plot_individual_spectra"] = plot_individual_spectra

        if linestyle is not self.add_entry.__defaults__[4]:  # if linestyle is provided
            df_temp["linestyle"] = linestyle

        if alpha is not self.add_entry.__defaults__[6]:  # if alpha is provided
            df_temp["alpha"] = alpha

        self.df = pd.concat([self.df, df_temp], ignore_index=True)
        self.df_fit_coeffs = pd.concat([self.df_fit_coeffs, df_temp_fit], ignore_index=True)

    def create_boxplot(self):
        # iterate over all group numbers
        medianprops = dict(linewidth=1.5, linestyle='-', color='k')
        sns.set_style("whitegrid")
        sns.set_palette("bright")
        #rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
        #rc('text', usetex=True)

        for i in self.df.group_number.drop_duplicates().tolist():
            df = self.df_fit_coeffs[self.df_fit_coeffs["group_number"] == i]
            # set the label: get the state and name from the df
            label = str(df["group_name"].iloc[0])

            fig, axes = plt.subplots(1, 4, figsize=(9.5 * 3 / 3, 4.2))
            fig.subplots_adjust(top=0.3)

            plt.subplot(131)
            label1 = str(np.round(np.mean(df["height"]), 3)) + " $\pm$ " + str(np.round(np.std(df["height"]), 3))

            plt.boxplot(df[["height"]].values, labels=[label1], medianprops=medianprops, showfliers=False)
            xs = np.random.normal(1, 0.04,
                                  df[["height"]].values.shape[0])  # adds jitter to the data points - can be adjusted
            plt.scatter(xs, df["height"], alpha=0.4, color="b")
            plt.ylabel("height [a.u.]")

            plt.subplot(132)
            label2 = str(np.round(np.mean(df["mean"]), 3)) + " $\pm$ " + str(
                np.round(np.std(df["mean"]), 3)) + " eV"

            plt.boxplot(df[["mean"]].values, labels=[label2], medianprops=medianprops, showfliers=False)
            xs = np.random.normal(1, 0.04,
                                  df[["mean"]].values.shape[0])  # adds jitter to the data points - can be adjusted
            plt.scatter(xs, df["mean"], alpha=0.4, color="r")
            plt.ylabel("center [eV]")

            plt.subplot(133)
            label3 = str(np.round(np.mean(df["fwhm"]), 3)) + " $\pm$ " + str(
                np.round(np.std(df["fwhm"]), 3)) + " eV"

            plt.boxplot(df[["fwhm"]].values, labels=[label3], medianprops=medianprops, showfliers=False)
            xs = np.random.normal(1, 0.04,
                                  df[["fwhm"]].values.shape[0])  # adds jitter to the data points - can be adjusted
            plt.scatter(xs, df["fwhm"], alpha=0.4, color="m")
            plt.ylabel("FWHM [eV]")

            plt.suptitle(
                label + " - " + str(np.shape(df)[0]) + " Particles")
            fig.tight_layout()
            fig.savefig(self.save_folder + self.timestr + "_Box_Plot_Group_" + str(i) + "_.pdf")
            fig.savefig(self.save_folder + self.timestr + "_Box_Plot_Group_" + str(i) + "_.png",dpi=600,
                        transparent=False)

        plt.close("all")

    def create_boxplot_side_by_side(self):
        # iterate over all group numbers
        medianprops = dict(linewidth=1.5, linestyle='-', color='k')
        sns.set_style("whitegrid")
        sns.set_palette("bright")
        #rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
        #rc('text', usetex=True)

        bplot_columns = ["height", "mean", "fwhm"]
        bplot_axname = ["height [a.u.]", "center [eV]", "FWHM [eV]"]
        bplot_units = ["[a.u]", "[eV]", "[eV]"]

        for column, axname, unit in zip(bplot_columns, bplot_axname, bplot_units):

            fig, axes = plt.subplots(1,  np.max(self.df.group_number), figsize=(9.5 / 3 * np.max(self.df.group_number), 4.2))

            ax = plt.gca()

            for i in self.df.group_number.drop_duplicates().tolist():
                df = self.df_fit_coeffs[self.df_fit_coeffs["group_number"] == i]
                # set the label: get the state and name from the df
                label = str(df["group_name"].iloc[0])

                fig.subplots_adjust(top=0.3)

                if i == 1:
                    ax = plt.subplot(1, np.max(self.df.group_number), i)
                else:
                    plt.subplot(1, np.max(self.df.group_number), i, sharey=ax)

                label2 = str(np.round(np.mean(df[column]), 3)) + " $\pm$ " + str(
                    np.round(np.std(df[column]), 3)) + " " + unit

                plt.boxplot(df[[column]].values, labels=[label2], medianprops=medianprops, showfliers=False)
                xs = np.random.normal(1, 0.04,
                                      df[[column]].values.shape[0])  # adds jitter to the data points - can be adjusted
                plt.scatter(xs, df[column], alpha=0.4, color="b")
                plt.ylabel(axname)
                plt.title(label + " No: " + str(np.shape(df)[0]))



            #plt.suptitle(
            #    label + " - " + str(np.shape(df)[0]) + " Particles - Group " + str(i))
            fig.tight_layout()
            fig.savefig(self.save_folder + self.timestr + "_" + column + "_Box_Plot_Side_by_Side.pdf")
            fig.savefig(self.save_folder + self.timestr + "_" + column + "_Box_Plot_Side_by_Side.png", dpi=600,
                        transparent=False)

        plt.close("all")

    def create_correlation_plots_side_by_side(self):
        # iterate over all group numbers
        medianprops = dict(linewidth=1.5, linestyle='-', color='k')
        sns.set_style("darkgrid")
        sns.set_palette("bright")
        #rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
        #rc('text', usetex=True)

        corr_columns_x = ["mean", "mean", "fwhm",
                          "particle_no", "particle_no", "particle_no"]
        bplot_axname_x = ["center [eV]", "center [eV]", "FWHM [eV]",
                          "measurement number", "measurement number", "measurement number"]
        corr_columns_y = ["height", "fwhm", "height",
                          "height", "mean", "fwhm"]
        bplot_axname_y = ["height [a.u.]", "FWHM [eV]", "height [a.u.]",
                          "height [a.u.]", "center [eV]", "FWHM [eV]"]

        for column_x, column_y, axname_x, axname_y in zip(corr_columns_x, corr_columns_y, bplot_axname_x, bplot_axname_y):

            fig, axes = plt.subplots(1,  np.max(self.df.group_number), figsize=(9.5 / 3 * np.max(self.df.group_number), 4.2))

            ax = plt.gca()

            for i in self.df.group_number.drop_duplicates().tolist():
                df = self.df_fit_coeffs[self.df_fit_coeffs["group_number"] == i]
                # set the label: get the state and name from the df
                label = str(df["group_name"].iloc[0])

                fig.subplots_adjust(top=0.3)

                if i == 1:
                    ax = plt.subplot(1, np.max(self.df.group_number), i)
                else:
                    # if you want shared x or y axes for these plots uncomment the following line
                    #plt.subplot(1, np.max(self.df.group_number), i, sharey=ax, sharex=ax)
                    plt.subplot(1, np.max(self.df.group_number), i)  # no shared axes

                #label2 = str(np.round(np.mean(df[column]), 3)) + " $\pm$ " + str(
                #    np.round(np.std(df[column]), 3)) + " eV"

                plt.plot(df[[column_x]].values, df[[column_y]].values, "o", color="b", alpha=0.3)
                plt.xlabel(axname_x)
                plt.ylabel(axname_y)

                plt.title(label + " No: " + str(np.shape(df)[0]))



            #plt.suptitle(
            #    label + " - " + str(np.shape(df)[0]) + " Particles - Group " + str(i))
            fig.tight_layout()
            fig.savefig(self.save_folder + self.timestr + "_" + column_x + "_" + column_y + "_Corr_Plot_Side_by_Side.pdf")
            fig.savefig(
                self.save_folder + self.timestr + "_" + column_x + "_" + column_y + "_Corr_Plot_Side_by_Side.png",
            dpi=600, transparent=False)

        plt.close("all")


    def export_all_csv(self):
        self.df.to_csv(self.save_folder + "csv/" + self.timestr + "_Averaged_Spectrum_Data.csv", index=False)
        self.df_fit_coeffs.to_csv(self.save_folder + "csv/" + self.timestr + "_Fit_Coeffs.csv", index=False)

    def average_spectra(self, coarse_grid=False, xlim=[]):
        """
        average the spectra by the previously (manually set) group number
        coarse_grid: If true, the x energies will be rounded first and then averaged
        """
        fig = plt.figure()  # create figure
        sns.set_style("darkgrid")
        sns.set_palette("bright")
        #rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
        #rc('text', usetex=True)
        y_lim_max = 0  # find the maximum y value in the used data
        print("0")

        # iterate over all group numbers
        for i in self.df.group_number.drop_duplicates().tolist():

            df_temp = self.df[self.df["group_number"] == i]

            if coarse_grid:  # apply rounding to energies, s.t. measurement from mutliple days can be combined
                df_temp["x"] = np.round(df_temp["x"].copy() / 4.2, 3) * 4.2

            df_group = df_temp.groupby("x", as_index=False)["y"].mean()  # group df here
            y_lim_max = y_lim_max if y_lim_max > np.max(df_group.y) else np.max(df_group.y)  # update ylim

            # set the label: get the name from the df
            label = str(df_temp["group_name"].iloc[0])

            alpha_fct_call = 0.1 if not ("alpha" in df_temp) or pd.isnull(df_temp["alpha"].iloc[0]) else df_temp["alpha"].iloc[0]
            if not(alpha_fct_call == 0.1): print("alpha provided")

            # plot the individual spectra
            if ("plot_individual_spectra" in df_temp) and not pd.isnull(df_temp["plot_individual_spectra"].iloc[0]):
                if ("color" in df_temp) and not pd.isnull(df_temp["color"].iloc[0]):
                    color = str(df_temp["color"].iloc[0])
                    for k in df_temp.counter_spectrum.drop_duplicates().tolist():
                        df_temp_spectrum = df_temp[df_temp["counter_spectrum"] == k]
                        plt.plot(df_temp_spectrum["x"], df_temp_spectrum["y"], ".", alpha=alpha_fct_call, color=color)

                else:
                    for k in df_temp.counter_spectrum.drop_duplicates().tolist():
                        df_temp_spectrum = df_temp[df_temp["counter_spectrum"] == k]
                        plt.plot(df_temp_spectrum["x"], df_temp_spectrum["y"], ".", alpha=alpha_fct_call)

            else: # average (do not plot the indiviudal spectra)
                if ("color" in df_temp) and not pd.isnull(df_temp["color"].iloc[0]) and \
                        ("linestyle" not in df_temp or (
                        pd.isnull(df_temp["linestyle"].iloc[0]))):  # color but no linestyle
                    print("1")
                    color = str(df_temp["color"].iloc[0])
                    plt.plot(df_group["x"], df_group["y"], "o", alpha=0.5, label=label, color=color)

                elif ("color" not in df_temp or pd.isnull(df_temp["color"].iloc[0])) and \
                        ("linestyle" in df_temp) and not pd.isnull(
                    df_temp["linestyle"].iloc[0]):  # only linestyle in short form
                    linestyle = str(df_temp["linestyle"].iloc[0])
                    print("2")
                    plt.plot(df_group["x"], df_group["y"], linestyle, alpha=0.5, label=label)

                elif ("color" in df_temp) and not pd.isnull(df_temp["color"].iloc[0]) and \
                        ("linestyle" in df_temp) and not pd.isnull(df_temp["linestyle"].iloc[0]):  # linestyle and color
                    linestyle = str(df_temp["linestyle"].iloc[0])
                    color = str(df_temp["color"].iloc[0])
                    print("3")
                    plt.plot(df_group["x"], df_group["y"], linestyle, alpha=0.5, label=label, color=color)

                else:
                    print("4")
                    plt.plot(df_group["x"], df_group["y"], "o", alpha=0.5, label=label)

        plt.xlabel("Emission Energy [eV]")
        plt.ylabel("counts [a.u.]")
        #plt.legend(loc=1)
        plt.legend()

        if xlim is not self.average_spectra.__defaults__[1]:  # xlim is provided
            plt.xlim(xlim)

        #plt.ylim([0, 1.15 * y_lim_max])
        fig.tight_layout()
        fig.savefig(self.save_folder + self.timestr + "_Averaged_Spectrum" + ".pdf")
        fig.savefig(self.save_folder + self.timestr + "_Averaged_Spectrum" + ".png", dpi=600, transparent=False)
        plt.close("all")

def create_folder(folderpath: str):
    """
    Create a folder given a path to a folder. This does not work recursively atm (i.e., you have to create every level of a
    nested folder step by step
    """
    if not (os.path.exists(folderpath)):  # if folder not found
        os.mkdir(folderpath)  # create path
        print("folder created: ", folderpath)
    return folderpath


if __name__ == "__main__":
    """
    Some of the add_entry options include:

    - "group_number": The group number is used to indicate which measurement days belong together; the spectra will be grouped by this 
        number. The numbers need to be positive integers starting at 1 in a connected intervall:
            Right: 1 2 3 4 ...
            Wrong: 0 1 2 3 ...
            Wrong: 1 3 4 5 ...
    - "group_name": String containing a user-chosen name for the specific group
    - "state": filter by state
        - state is in ["singular", "non_singular", "discarded"] these are based on g2 if available    
        - state is in [True, False] (careful: boolean values!) these are based on the calculations in Fitting_SingleDots
            and are therefore based on the threshold for chi-squared in that script (most likely: 1)
    - "threshold_chisqr" : filter by chi-squared value: this implements (almost) excactly the same filtering as in the Fitting_SingleDots 
            scripts and can therefore be used to quickly test different thresholds! Values need to be larger than 1
            with sensible values in the range 0.1 ... 1.0.
    - "selected_blobs" select specfic spectra: here an array can be entered with the Andor File numbers of spectra one 
        wants to visualize
        
        
    Plotting Options (these might not always work 100% correctly). These options are mainly for the spectra plots
    - "color": specifiy any color. The shorthand notation for the colors  ("c", "m", "y", "k", "r", "g", "b") is mapped
        to match the specified seaborn color palette. It is also possible to provide a hex value. Refer to the official
        documentation
    - "linestyle": different markers can be selected ("o", "s", "<", ".")
    
    ATTENTION: Currently only the lorentzian fit data is used! #TODO
    """

    save_folder_path_global = create_folder("T:/Gabriel/Data/23-01-26_AP_3_120_MultipleSingleDot/andor_export/comparison/")
    _ = create_folder(save_folder_path_global + "csv/")

    Import = csvFiles()

    Import.add_entry(
        folder="T:/Gabriel/Data/23-01-26_AP_3_120_MultipleSingleDot/andor_export/",
        group_number=1, group_name="AP-3-120", color="b", threshold_chisqr=1, upper_limit_fwhm=0.1)

    Import.create_boxplot()
    Import.create_boxplot_side_by_side()
    Import.create_correlation_plots_side_by_side()

    Import.average_spectra(coarse_grid=False)
    Import.export_all_csv()
