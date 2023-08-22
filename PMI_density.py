# Hillel Darshan 8/22/23
import numpy as np
import pandas as pd
import plotly.graph_objects as go

























if __name__ == '__main__':
    files_location = "../PMI_analysis/working_data_PMI/PMI_Adi_DATA/"
    df = pd.read_csv(files_location + "RUSH_patients_info.csv")
    draw_volcano_plots(df)

    # check spearman correlation
    # spearman_rin_pmi(df)

    print("All Done")