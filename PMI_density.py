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


def calculate_patients_df():
    # df of patients
    pmi_projid = pd.read_csv(working_dir + "rush dataset_662_basic_08-10-2023.csv")
    cold_projid = pd.read_csv(working_dir + "cold_data_ros.csv")
    lab_id = pd.read_csv(working_dir + "lab_name_id_table.csv")
    lab_id = lab_id.merge(pmi_projid, on='projid', how='left')
    lab_id = lab_id.merge(cold_projid, on='projid', how='left')

    # keep only the cols we want
    cols_to_drop = list(lab_id.columns)
    cols_to_drop = [e for e in cols_to_drop if e not in ('Name', 'projid', 'pmi', 'cogdx', 'msex',
                                                         'ad_reagan', 'braaksc')]
    lab_id.drop(cols_to_drop, axis=1, inplace=True)

    return lab_id


def calculate_trf_df():
    # df of tRFs
    trf_df = pd.read_csv(working_dir + "tRNA_Exclusive_Combined_data.csv")
    trf_df.columns = ['tRF name'] + list(trf_df.columns[1:])
    meta = pd.read_csv(working_dir + "tRF_meta.csv")
    meta.columns = ['tRF name'] + list(meta.columns[1:])
    trf_df = trf_df.merge(meta, on='tRF name', how='left')
    trf_df.drop(['trf', 'tRF sequence', 'details', 'trna', 'nuclear'], axis=1, inplace=True)

    return trf_df


if __name__ == '__main__':
    working_dir = "../PMI_analysis/working_data_PMI/PMI_Adi_DATA/"
    df = calculate_trf_df()
    patients_df = calculate_patients_df()
    sums_df = calculate_sums_df(df)

    sums_df_t = sums_df.T
    new_columns = [f"{sums_df_t.iloc[0, i]}_{sums_df_t.iloc[1, i]}" for i in range(len(sums_df_t.columns))]
    sums_df_t.columns = new_columns
    sums_df_t = sums_df_t.iloc[2:]
    sums_df_t['Name'] = sums_df_t.index
    sums_df_t['Name'] = sums_df_t['Name'].str.split('_').str[:2].str.join(' ')
    sums_df_t = sums_df_t.merge(patients_df, on='Name', how='left', sort=False)

    intercept_lst = [
                     ('AD pmi <=5', (sums_df_t['cogdx'] == 4) & (sums_df_t['pmi'] <= 5)),
                     ('AD 5<pmi<=7', (sums_df_t['cogdx'] == 4) & (sums_df_t['pmi'] <= 7)  & (sums_df_t['pmi'] > 5)),
                     ('AD 7<pmi<=9', (sums_df_t['cogdx'] == 4) & (sums_df_t['pmi'] <= 9)  & (sums_df_t['pmi'] > 7)),
                     ('AD 9<pmi<=11', (sums_df_t['cogdx'] == 4) & (sums_df_t['pmi'] <= 11)  & (sums_df_t['pmi'] > 9)),
                     ('AD 11<pmi<=13', (sums_df_t['cogdx'] == 4) & (sums_df_t['pmi'] <= 13)  & (sums_df_t['pmi'] > 11)),
                     ('AD 13<pmi<=15', (sums_df_t['cogdx'] == 4) & (sums_df_t['pmi'] <= 15)  & (sums_df_t['pmi'] > 13)),
                     ('AD 15<pmi', (sums_df_t['cogdx'] == 4) & (sums_df_t['pmi'] > 15)),

                     ('MMCI pmi <=5', (sums_df_t['cogdx'] == 2) & (sums_df_t['pmi'] <= 5)),
                     ('MMCI 5<pmi<=7', (sums_df_t['cogdx'] == 2) & (sums_df_t['pmi'] <= 7) & (sums_df_t['pmi'] > 5)),
                     ('MMCI 7<pmi<=9', (sums_df_t['cogdx'] == 2) & (sums_df_t['pmi'] <= 9) & (sums_df_t['pmi'] > 7)),
                     ('MMCI 9<pmi<=11', (sums_df_t['cogdx'] == 2) & (sums_df_t['pmi'] <= 11) & (sums_df_t['pmi'] > 9)),
                     ('MMCI 11<pmi<=13', (sums_df_t['cogdx'] == 2) & (sums_df_t['pmi'] <= 13) & (sums_df_t['pmi'] > 11)),
                     ('MMCI 13<pmi<=15', (sums_df_t['cogdx'] == 2) & (sums_df_t['pmi'] <= 15) & (sums_df_t['pmi'] > 13)),
                     ('MMCI 15<pmi', (sums_df_t['cogdx'] == 2) & (sums_df_t['pmi'] > 15)),

                     ('NCI pmi <=5', (sums_df_t['cogdx'] == 1) & (sums_df_t['pmi'] <= 5)),
                     ('NCI 5<pmi<=7', (sums_df_t['cogdx'] == 1) & (sums_df_t['pmi'] <= 7) & (sums_df_t['pmi'] > 5)),
                     ('NCI 7<pmi<=9', (sums_df_t['cogdx'] == 1) & (sums_df_t['pmi'] <= 9) & (sums_df_t['pmi'] > 7)),
                     ('NCI 9<pmi<=11', (sums_df_t['cogdx'] == 1) & (sums_df_t['pmi'] <= 11) & (sums_df_t['pmi'] > 9)),
                     ('NCI 11<pmi<=13', (sums_df_t['cogdx'] == 1) & (sums_df_t['pmi'] <= 13) & (sums_df_t['pmi'] > 11)),
                     ('NCI 13<pmi<=15', (sums_df_t['cogdx'] == 1) & (sums_df_t['pmi'] <= 15) & (sums_df_t['pmi'] > 13)),
                     ('NCI 15<pmi', (sums_df_t['cogdx'] == 1) & (sums_df_t['pmi'] > 15))
                     ]
    for intercept in intercept_lst:
        new_row = sums_df_t[intercept[1]].sum()
        new_row.name = intercept[0]
        sums_df_t = sums_df_t.append(new_row)
    print("All Done")