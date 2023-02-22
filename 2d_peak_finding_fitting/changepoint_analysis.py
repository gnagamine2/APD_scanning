import pandas as pd
import numpy as np
import ruptures as rpt
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import statistics
import os
import re
import math
sns.set()
sns.set_palette("bright")
palette = sns.color_palette("bright")
matplotlib.use("Qt5Agg")

# These are without the iterations for singular arrays
"""
def partition_js(data):
    data = sorted(data)

    best_split_indices = [len(data)]
    best_split_score = score(best_split_indices, data)

    for split_index in range(1, len(data) - 2):
        split_indices = [split_index, len(data)]
        s = score(split_indices, data)
        # print("s overall", s)
        if s > best_split_score:
            best_split_score = s
            best_split_indices = split_indices

    for first_split in range(1, len(data) - 3):
        for second_split in range(first_split + 2, len(data) - 2):
            split_indices = [first_split, second_split, len(data)]
            # print(split_indices)
            s = score(split_indices, data)
            if s > best_split_score:
                best_split_score = s
                best_split_indices = split_indices

    return len(best_split_indices), best_split_indices, best_split_score


def score(splits, data: list[float]) -> float:
    score = 0
    split_start = 0
    for i in splits:
        # print("i", i)
        group = data[split_start: i + 1]
        # print("group", group)
        split_start = i + 1
        # print("variance: ", statistics.variance(group))
        score += statistics.variance(group)
        # instead of variance median absolute distance to median
        # https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
        # print("score", score)
    return -(len(splits)**(2.8) + 1) * score
"""

def partition_js(data):
    data = sorted(data)

    best_split_indices = [len(data)]
    best_split_score = score(best_split_indices, data)
    #print(" ")

    for split_index in range(0, len(data) - 1):
        split_indices = [split_index, len(data)]
        s = score(split_indices, data)
        # print("s overall", s)
        if s > best_split_score:
            best_split_score = s
            best_split_indices = split_indices
        #print(" ")

    for first_split in range(0, len(data) - 1):
        for second_split in range(first_split + 1, len(data) - 1):
            split_indices = [first_split, second_split, len(data)]
            # print(split_indices)
            s = score(split_indices, data)
            if s > best_split_score:
                best_split_score = s
                best_split_indices = split_indices
            #print(" ")
        #print(" ")

    return len(best_split_indices), best_split_indices, best_split_score

def score(splits, data: list[float]) -> float:
    score = 0
    split_start = 0
    for i in splits:
        # print("i", i)
        group = data[split_start: i + 1]
        #print("group", group)
        split_start = i + 1
        # print("variance: ", statistics.variance(group))
        if len(group) == 1:
            score += 0
        else:
            score += statistics.variance(group)
        # instead of variance median absolute distance to median
        # https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
    #print("score", score)
    return -(len(splits) + 1) * score

def changepoint_analysis(filename: str):
    # ATTENTION: USE EVERY THIRD POINTS
    df = pd.DataFrame()
    df["x"] = df_temp.x[::3]
    df["y"] = df_temp.y[::3]

    points = df.y.values

    # RUPTURES PACKAGE
    # Changepoint detection with the Pelt search method
    model = "rbf"
    algo = rpt.Pelt(model=model).fit(points)
    result = algo.predict(pen=1)
    # print("result")

    # pad result with the start value of df.x
    result_array = np.insert(result, 0, np.min(df.x), axis=0)

    count_levels = np.zeros(np.size(result_array) - 1)
    time_count = np.zeros(np.size(result_array) - 1)
    # df_res = pd.DataFrame()
    # TODO: Check if lowest count level is at bg
    # TODO: Maybe put a weight on the count states? I.e., longer more important

    for i in range(0, np.size(result_array) - 1):
        i_start = result_array[i]
        i_end = result_array[i + 1] - 1
        count_levels[i] = np.median(df.y[i_start:i_end])
        time_count[i] = i_end - i_start
        # length_of_interval.append(i_end - i_start)
        # df_res.at[i, 'median_levels'] = np.median(df.y[i_start:i_end])
        # df_res.at[i, 'i_start'] = i_start
        # df_res.at[i, 'i_end'] = i_end
        # df_res.at[i, 't_start'] = df.x[i_start]
        # df_res.at[i, 't_end'] = df.x[i_end]
        # df_res.at[i, 'time_diff'] = df.x[i_end] - df.x[i_start]
    # print("iterated")
    # print(df_res.median_levels)
    time_count = time_count * np.max(df.x) / np.max(result_array)

    time_below_threshold = 0
    for (time, level) in zip(time_count, count_levels):
        if level < 100:
            time_below_threshold += time
    time_below_threshold

    df_histogram.at[particle_no, "time_below_threshold"] = time_below_threshold

    if len(count_levels) == 1:
        array.append((particle_no, 1))
        if (df_g2[df_g2["particle_no"] == particle_no].state == "singular").values[0]:
            array_correct.append(False)
            array_correct_no.append((particle_no, 1))

    no_states, _, _ = partition_js(count_levels.tolist())
    array.append((particle_no, no_states))
    # print("partition completed")

    if (df_g2[df_g2["particle_no"] == particle_no].state == "singular").values[0]:
        if no_states == 2:
            array_correct.append(True)
            array_correct_no.append((particle_no, no_states))
        else:
            array_correct.append(False)
            array_wrong_no.append((particle_no, no_states))
    elif no_states == 2:
        array_wrong_no.append((particle_no, no_states))

def histogram_integral_max(df: pd.DataFrame(), bright_dark_cutoff: int, filename: str):
    """
    Calculate the integral based on a count threshold for the dark and bright states. Find the maximum within these
    regions.
    """
    # create the histogram
    histogram = np.histogram(df.y.values, int(max(df.y.values / 2)), [0, int(max(df.y.values))])
    histogram_x = histogram[0]
    histogram_y = 0.5 * (histogram[1][:-1] + histogram[1][1:])

    # check if counts includes the bright_dark_cutoff value
    bright_dark_cutoff_idx = (
        np.abs(histogram_y - bright_dark_cutoff)).argmin()  # find the idx closest to bright_dark_cutoff

    dark_sum = np.sum(histogram_x[:bright_dark_cutoff_idx])
    bright_sum = np.sum(histogram_x[bright_dark_cutoff_idx:])

    dark_max = np.max(histogram_x[:bright_dark_cutoff_idx])
    bright_max = np.max(histogram_x[bright_dark_cutoff_idx:])

    plt.figure()
    plt.plot(histogram_x[:bright_dark_cutoff_idx], histogram_y[:bright_dark_cutoff_idx])
    plt.plot(histogram_x[bright_dark_cutoff_idx:], histogram_y[bright_dark_cutoff_idx:])
    plt.xlabel('occurrence')
    plt.ylabel('counts / %i ms' % (0.1 * 1e3))  # TODO: Careful! This bin value is fixed
    plt.savefig("Data/Histogram/" + filename[:-4] + "_Histogram.pdf")

    # add values to global (!) df
    df_histogram.at[particle_no, "dark_sum"] = dark_sum
    df_histogram.at[particle_no, "bright_sum"] = bright_sum
    df_histogram.at[particle_no, "dark_max"] = dark_max
    df_histogram.at[particle_no, "bright_max"] = bright_max

    return dark_sum, bright_sum, dark_max, bright_max

if __name__ == "__main__":
    df_g2 = pd.read_csv("Data/state_files.csv")
    df_g2['particle_no'] = df_g2['name'].apply(
        lambda f: int(re.search(r"Particle_(\d{1,3})_", f).group(1)))
    lst_singular = df_g2[df_g2.state == "singular"].particle_no.tolist()
    lst_non_singular = df_g2[df_g2.state == "non_singular"].particle_no.tolist()



    count = 0
    array = []
    array_correct = []
    array_correct_no = [] #contains tuple of files identified correctly
    array_wrong_no = []
    df_histogram = pd.DataFrame()
    for filename in sorted(os.listdir("data/")):
        if not filename.endswith("csvIntTrace.csv"):
            continue

        print(filename)

        particle_no = int(re.search(r"Particle_(\d{1,3})_", filename).group(1))
        print("count: ", count)
        df_temp = pd.read_csv("data/" + filename)

        # determine the g2 value at zero delay
        g2_tau_zero = df_g2[df_g2["particle_no"] == particle_no]["g2_tau_zero"].values[0]
        # add g2 at zero to df
        df_histogram.at[particle_no, "g2_tau_zero"] = g2_tau_zero

        changepoint_analysis(filename);

        #histogram_integral_max(df=df_temp, bright_dark_cutoff=100, filename=filename)

        count += 1

    """
    # Create Histogram Plot
    df_temp = df_histogram
    df_temp.dropna(inplace=True)
    df_temp["sum_ratio"] = df_temp["bright_sum"] / df_temp["dark_sum"]
    df_temp["max_ratio"] = df_temp["bright_max"] / df_temp["dark_max"]

    df_sum = df_temp[df_temp["sum_ratio"] < 100] #outlier rejection

    fig = plt.figure()
    sns.scatterplot(data=df_sum, x="g2_tau_zero", y="sum_ratio", label="all")
    sns.scatterplot(data=df_sum[df_sum["sum_ratio"] < 7], x="g2_tau_zero", y="sum_ratio", label="sum_ratio < 7")
    fig.savefig("Data/g2_vs_sumRatio.pdf")

    fig = plt.figure()
    sns.scatterplot(data=df_temp, x="g2_tau_zero", y="max_ratio", label="all")
    sns.scatterplot(data=df_sum[df_sum["max_ratio"] < 2], x="g2_tau_zero", y="max_ratio", label="max_ratio < 2")
    fig.savefig("Data/g2_vs_maxRatio.pdf")
    """


