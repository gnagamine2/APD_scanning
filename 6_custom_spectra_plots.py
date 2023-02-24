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

    def add_entry(self, filepath, color=None):
        """
        Loads data from a csv file and adds it as a new entry to the csvFiles object.
    
        Parameters
        ----------
        filepath : str
            Path to the csv file containing the data.
        color : str, optional
            The color to be used for this entry in the boxplot graph. If not provided,
            a default color will be used.
        """
        if color is None:
            color = self.default_color
    
        data = np.genfromtxt(filepath, delimiter=",")
        self.samples.append(data)
        self.colors.append(color)

    

    def create_boxplot(self, show_outliers=False, show_means=False, color_by_sample=False):
        """
        Creates a boxplot graph from the data in the csvFiles object.
    
        Parameters
        ----------
        show_outliers : bool, optional
            Whether or not to show outliers in the graph. Default is False.
        show_means : bool, optional
            Whether or not to show means in the graph. Default is False.
        color_by_sample : bool, optional
            Whether to color each box by sample (True) or to use a single color for all
            boxes (False). Default is False.
        """
        data = self.samples
        colors = self.colors
        num_samples = len(data)
    
        # Determine the box colors
        if color_by_sample:
            box_colors = colors
        else:
            box_colors = [self.default_color] * num_samples
    
        # Create the boxplot graph
        fig, ax = plt.subplots()
        ax.boxplot(data, showfliers=show_outliers, showmeans=show_means,
                   boxprops={'linewidth': 2, 'color': 'k'},
                   whiskerprops={'linewidth': 2, 'color': 'k'},
                   capprops={'linewidth': 2, 'color': 'k'},
                   medianprops={'linewidth': 2, 'color': 'r'},
                   meanprops={'linewidth': 2, 'color': 'b'},
                   patch_artist=True,
                   box_colors=box_colors)
    
        # Add labels and legend
        ax.set_xticklabels([f"Sample {i+1}" for i in range(num_samples)])
        ax.set_ylabel("Intensity (a.u.)")
        if color_by_sample:
            legend_handles = [mpatches.Patch(color=color, label=f"Sample {i+1}") for i, color in enumerate(colors)]
            ax.legend(handles=legend_handles, loc="upper right")
    
        # Save the graph
        if not os.path.exists(self.plot_folder):
            os.makedirs(self.plot_folder)
        fig.savefig(os.path.join(self.plot_folder, "boxplot.pdf"), bbox_inches="tight")

        
        

    # def create_boxplot(self):
    #     # iterate over all group numbers
    #     medianprops = dict(linewidth=1.5, linestyle='-', color='k')
    #     sns.set_style("whitegrid")
    #     sns.set_palette("bright")
    #     #rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    #     #rc('text', usetex=True)

    #     for i in self.df.group_number.drop_duplicates().tolist():
    #         df = self.df_fit_coeffs[self.df_fit_coeffs["group_number"] == i]
    #         # set the label: get the state and name from the df
    #         label = str(df["group_name"].iloc[0])

    #         fig, axes = plt.subplots(1, 4, figsize=(9.5 * 3 / 3, 4.2))
    #         fig.subplots_adjust(top=0.3)

    #         plt.subplot(131)
    #         label1 = str(np.round(np.mean(df["height"]), 3)) + " $\pm$ " + str(np.round(np.std(df["height"]), 3))

    #         plt.boxplot(df[["height"]].values, labels=[label1], medianprops=medianprops, showfliers=False)
    #         xs = np.random.normal(1, 0.04,
    #                               df[["height"]].values.shape[0])  # adds jitter to the data points - can be adjusted
    #         plt.scatter(xs, df["height"], alpha=0.4, color="b")
    #         plt.ylabel("height [a.u.]")

    #         plt.subplot(132)
    #         label2 = str(np.round(np.mean(df["mean"]), 3)) + " $\pm$ " + str(
    #             np.round(np.std(df["mean"]), 3)) + " eV"

    #         plt.boxplot(df[["mean"]].values, labels=[label2], medianprops=medianprops, showfliers=False)
    #         xs = np.random.normal(1, 0.04,
    #                               df[["mean"]].values.shape[0])  # adds jitter to the data points - can be adjusted
    #         plt.scatter(xs, df["mean"], alpha=0.4, color="r")
    #         plt.ylabel("center [eV]")

    #         plt.subplot(133)
    #         label3 = str(np.round(np.mean(df["fwhm"]), 3)) + " $\pm$ " + str(
    #             np.round(np.std(df["fwhm"]), 3)) + " eV"

    #         plt.boxplot(df[["fwhm"]].values, labels=[label3], medianprops=medianprops, showfliers=False)
    #         xs = np.random.normal(1, 0.04,
    #                               df[["fwhm"]].values.shape[0])  # adds jitter to the data points - can be adjusted
    #         plt.scatter(xs, df["fwhm"], alpha=0.4, color="m")
    #         plt.ylabel("FWHM [eV]")

    #         plt.suptitle(
    #             label + " - " + str(np.shape(df)[0]) + " Particles")
    #         fig.tight_layout()
    #         fig.savefig(self.save_folder + self.timestr + "_Box_Plot_Group_" + str(i) + "_.pdf")
    #         fig.savefig(self.save_folder + self.timestr + "_Box_Plot_Group_" + str(i) + "_.png",dpi=600,
    #                     transparent=False)

    #     plt.close("all")

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

    save_folder_path_global = create_folder("C:/DataAnalysis/AP_3_108/plots/")
    _ = create_folder(save_folder_path_global + "csv/")

    Import = csvFiles()

    Import.add_entry(
        folder="C:/DataAnalysis/AP_3_108/22-11-08_AP-3-108_MultipleDots_2Batches/andor_export/",
        group_number=1, group_name="AP-3-108", color="b", threshold_chisqr=1, upper_limit_fwhm=0.15)
    Import.add_entry(
        folder="C:/DataAnalysis/AP_3_108/22-11-11-AP-3-108_MultipleSingleDot_Measurement/andor_export/",
        group_number=1, group_name="AP-3-108", color="b", threshold_chisqr=1, upper_limit_fwhm=0.15)
    Import.add_entry(
        folder="C:/DataAnalysis/AP_3_108/22-11-18-AP-3-108_MultipleSingleDot_Measurement/andor_export/",
        group_number=1, group_name="AP-3-108", color="b", threshold_chisqr=1, upper_limit_fwhm=0.15)
    Import.add_entry(
        folder="C:/DataAnalysis/AP_3_108/22-11-23-AP_3_108/andor_export_2/",
        group_number=1, group_name="AP-3-108", color="b", threshold_chisqr=1, upper_limit_fwhm=0.15)


    color_list = ['red', 'green', 'orange', 'purple']
    Import.add_entry("csv/3660-3700.csv", color_list=color_list)
    Import.add_entry("csv/3725-3765.csv", color_list=color_list)
    Import.add_entry("csv/3790-3830.csv", color_list=color_list)
    Import.add_entry("csv/3860-3900.csv", color_list=color_list)

    Import.create_boxplot(show_outliers=True, show_means=True, color_by_sample=True)



    # Import.create_boxplot()
    # Import.create_boxplot_side_by_side()
    # Import.create_correlation_plots_side_by_side()

    # Import.average_spectra(coarse_grid=False)
    # Import.export_all_csv()
