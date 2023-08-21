# Hillel Darshan 8/21/23
import numpy as np
import pandas as pd
import plotly.graph_objects as go

def draw_volcano_plots(df):
    # df = df[df['brain_region'] != "NucAcc"]
    df['brain_region_color'] = df['brain_region'].apply(lambda x: "red" if x == "NucAcc" else "blue")

    pl = go.Figure(
        go.Scatter(x=df["pmi"],
                   y=df["RIN"],
                   mode="markers",
                   marker=dict(color=df['brain_region_color'])))
    pl.update_layout(title="RIN / PMI by brain region",
                     xaxis_title="PMI",
                     yaxis_title="RIN")

    pl.show()




if __name__ == '__main__':
    files_location = "../PMI_analysis/working_data_PMI/PMI_Dana_DATA/"
    df = pd.read_csv(files_location + "RUSH_patients_info.csv")
    # draw_volcano_plots(df)

    # check spearman correlation
    df = df[df['brain_region'] != "NucAcc"]
    cols_to_keep = ['pmi', 'RIN']
    sp = df[cols_to_keep]
    sp = sp.corr(method='spearman', min_periods=1, numeric_only=False)


    print("All Done")
