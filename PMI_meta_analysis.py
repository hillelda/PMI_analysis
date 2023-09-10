# Hillel Darshan 8/21/23
# meta analysis, for understanding the features and their correlation to RIN and PMI

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.stats import spearmanr

def draw_volcano_plots(df1):
    # df1 = df1[df1['brain_region'] != "NucAcc"]
    df1 = df1[df1['study'] != "MAP"]
    df1['brain_region_color'] = df1['brain_region'].apply(lambda x: "red" if x == "NucAcc" else "blue")

    pl = go.Figure(
        go.Scatter(x=df1["pmi"],
                   y=df1["RIN"],
                   mode="markers",
                   marker=dict(color=df1['brain_region_color'])))
    pl.update_layout(title="RIN / PMI by brain region",
                     xaxis_title="PMI",
                     yaxis_title="RIN")

    pl.show()


def spearman_rin_pmi(data_frame):
    # data_frame = data_frame[data_frame['brain_region'] != "NucAcc"]
    # data_frame = data_frame[data_frame['msex'] != "Male"]
    data_frame = data_frame[data_frame['study'] == "MAP"]
    cols_to_keep = ['pmi', 'RIN']
    sp = data_frame[cols_to_keep]
    # sp1 = sp.corr(method='spearman', min_periods=1, numeric_only=False)
    sp2 = spearmanr(sp['pmi'], sp['RIN'], axis=1)
    print(sp2)


if __name__ == '__main__':
    files_location = "../PMI_analysis/working_data_PMI/PMI_Dana_DATA/"
    df = pd.read_csv(files_location + "RUSH_patients_info.csv")
    # raw_counts = pd.read_csv(files_location + "raw_counts_combined_tRF.csv")
    draw_volcano_plots(df)

    # check spearman correlation
    # spearman_rin_pmi(df)

    print("All Done")
