# Hillel Darshan 8/22/23
from itertools import product

import numpy as np
import pandas as pd
import plotly.graph_objects as go


def calculate_sums_df(patient_df):
    # Generate all possible combinations of conditions
    tRF_groups = patient_df['tRF type(s)'].unique()
    lens = patient_df['len'].unique()
    condition_combinations = list(product(tRF_groups, lens))
    patient_columns = list(patient_df.columns)[2:-2]
    # calculate sums
    cur_sums_df = pd.DataFrame(columns=['tRF type(s)', 'len'] + patient_columns)
    for condition_combination in condition_combinations:
        trf_type, trf_len = condition_combination
        sum_condition = (patient_df['tRF type(s)'] == trf_type) & (patient_df['len'] == trf_len)

        if sum_condition.any():  # Check if there are matching rows
            sum_result = patient_df[sum_condition][patient_columns].sum()
            cur_sums_df.loc[len(cur_sums_df)] = [trf_type, trf_len] + list(sum_result)

    return cur_sums_df

if __name__ == '__main__':
    working_dir = "../PMI_analysis/working_data_PMI/PMI_Adi_DATA/"

    # df of tRFs
    df = pd.read_csv(working_dir + "tRNA_Exclusive_Combined_data.csv")
    df.columns = ['tRF name'] + list(df.columns[1:])
    meta = pd.read_csv(working_dir + "tRF_meta.csv")
    meta.columns = ['tRF name'] + list(meta.columns[1:])
    df = df.merge(meta, on='tRF name', how='left')
    df.drop(['trf', 'tRF sequence', 'details', 'trna', 'nuclear'] , axis=1, inplace=True)

    # df of patients
    pmi_projid = pd.read_csv(working_dir + "rush dataset_662_basic_08-10-2023.csv")
    lab_id = pd.read_csv(working_dir + "lab_name_id_table.csv")
    lab_id = lab_id.merge(pmi_projid, on='projid', how='left')
    lab_id.drop(['apoe_genotype', 'study'], axis=1, inplace=True)

    sums_df = calculate_sums_df(df)

    print("All Done")