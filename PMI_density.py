# Hillel Darshan 8/22/23

from itertools import product
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ########### Global vars ########
WORKING_DIR = "../PMI_analysis/working_data_PMI/PMI_Adi_DATA/"
DATA_SET = "rush dataset_662_basic_08-10-2023.csv"
COLD_DATA = "cold_data_ros.csv"
LAB_ID_FILE = "lab_name_id_table.csv"
TRF_FILE = "tRNA_Exclusive_Combined_data.csv"
TRF_META = "tRF_meta.csv"

PNG_OUTPUT = '../PMI_analysis/density_python.png'
OUTPUT_FILE = '../PMI_analysis/res_summed_by_type_len.csv'
# ################################

def calculate_sums_df(patient_df):
    # Generate all possible combinations of conditions
    tRF_groups = patient_df['tRF type(s)'].unique()
    lens = patient_df['len'].unique()
    condition_combinations = list(product(tRF_groups, lens))
    patient_columns = list(patient_df.columns)[2:-2]
    # calculate sums
    sums_df = pd.DataFrame(columns=['tRF type(s)', 'len'] + patient_columns)
    for condition_combination in condition_combinations:
        trf_type, trf_len = condition_combination
        sum_condition = (patient_df['tRF type(s)'] == trf_type) & (patient_df['len'] == trf_len)

        if sum_condition.any():  # Check if there are matching rows
            sum_result = patient_df[sum_condition][patient_columns].sum()
            sums_df.loc[len(sums_df)] = [trf_type, trf_len] + list(sum_result)

    return sums_df


def calculate_patients_df():
    # df of patients
    pmi_projid = pd.read_csv(WORKING_DIR + DATA_SET)
    cold_projid = pd.read_csv(WORKING_DIR + COLD_DATA)
    lab_id = pd.read_csv(WORKING_DIR + LAB_ID_FILE)
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
    trf_df = pd.read_csv(WORKING_DIR + TRF_FILE)
    trf_df.columns = ['tRF name'] + list(trf_df.columns[1:])
    meta = pd.read_csv(WORKING_DIR + TRF_META)
    meta.columns = ['tRF name'] + list(meta.columns[1:])
    trf_df = trf_df.merge(meta, on='tRF name', how='left')
    trf_df.drop(['trf', 'tRF sequence', 'details', 'trna', 'nuclear'], axis=1, inplace=True)

    return trf_df


def calculate_result_df(sums_df):
    # transpose matrix, and re-name cols to be names of patients
    sums_df_t = sums_df.T
    new_columns = [f"{sums_df_t.iloc[0, i]}_{sums_df_t.iloc[1, i]}" for i in range(len(sums_df_t.columns))]
    sums_df_t.columns = new_columns
    sums_df_t = sums_df_t.iloc[2:]
    sums_df_t['Name'] = sums_df_t.index.str.split('_').str[:2].str.join(' ')
    sums_df_t = sums_df_t.merge(patients_df, on='Name', how='left', sort=False)

    # generate intercept list
    intercept_lst = []

    categories = [
        ('AD', 4),
        ('MMCI', 2),
        ('NCI', 1)
    ]

    pmi_ranges = [
        (0, 5),
        (5, 7),
        (7, 9),
        (9, 11),
        (11, 13),
        (13, 15),
        (15, 100)
    ]

    for category, cogdx in categories:
        for pmi_range in pmi_ranges:
            label = f'{category} {pmi_range[0]}<pmi<{pmi_range[1]}'
            condition = (sums_df_t['cogdx'] == cogdx) & (sums_df_t['pmi'] <= pmi_range[1])

            if pmi_range[0] > 0:
                condition &= (sums_df_t['pmi'] > pmi_range[0])

            intercept_lst.append((label, condition))

    res_dfs = []
    for intercept in intercept_lst:
        new_row = sums_df_t[intercept[1]].sum()
        new_row.name = intercept[0]
        res_dfs.append(new_row)
    result_df = pd.concat(res_dfs, axis=1)

    result_df = result_df.iloc[:-7] # delete last rows, sums of meta-data
    result_df['type'] = result_df.index.str.split('_').str[:1].str.join(' ')
    result_df['len'] = result_df.index.str.split('_').str[1]

    return result_df


def plot_data(df):
    # Normalize columns
    for col in range(0, df.shape[1] - 2):
        df.iloc[:, col] = df.iloc[:, col] / (df.iloc[:, col].sum() + 1)
    # Melt the dataframe
    res = pd.melt(df, id_vars=['type', 'len'], var_name='group')
    # Convert 'len' to numeric, handling any non-numeric values
    res['len'] = pd.to_numeric(res['len'], errors='coerce')
    res = res.dropna(subset=['len'])
    # Extract cond and pmi
    res['condition'] = res['group'].apply(lambda x: x.split(' ')[0])
    res['pmi'] = res['group'].apply(lambda x: x.split('<')[-1])
    # plot
    g = sns.FacetGrid(res, col='condition', row='type', margin_titles=True,
                      xlim=(res['len'].min(), res['len'].max()), ylim=(res['value'].min(), res['value'].max()))
    g.map_dataframe(sns.lineplot, x='len', y='value', hue='pmi', linewidth=0.5)
    g.add_legend(title='pmi', label_order=res['pmi'].unique().tolist())

    return g


if __name__ == '__main__':
    df = calculate_trf_df()
    patients_df = calculate_patients_df()
    cur_sums_df = calculate_sums_df(df)
    cur_res_df = calculate_result_df(cur_sums_df)
    cur_res_df.to_csv(OUTPUT_FILE)
    p = plot_data(cur_res_df)
    plt.savefig(PNG_OUTPUT, dpi=900)
    plt.show()
    print("All Done")