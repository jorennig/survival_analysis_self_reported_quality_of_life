import pandas as pd

def time_to_event(df):
    df['shift_1'] = df['event'].shift(-1)
    df['shift_2'] = df['event'].shift(-2)
    
    df['change'] = (df['event']==df['shift_1']) & \
                   (df['event']==df['shift_2'])
    
    if df['change'].eq(False).all():
        weeks_to_event = df['weeks'].max()
        event_detected = False
    else:
        idx = df['change'][df['change']==True].index[0]
        weeks_to_event = df.loc[idx,]['weeks']
        event_detected = True
    return pd.Series([weeks_to_event, event_detected], index=['weeks_to_event', 'event_detected'])

input_file = snakemake.input.data
mdc_file = snakemake.input.mdc
output_data = snakemake.output[0]
features = snakemake.params.features
# input_file = r'\\exports.hps.kau.science.roche.com\pred\rpmda\BN40423-v2\jira\RPMDA-7138-longitudinal-analysis\data\features\STANDARD\0\demographics.tsv'
# mdc_file = r'\\exports.hps.kau.science.roche.com\pred\rpmda\BN40423-v2\jira\RPMDA-7438-cross-study-median-mdc-values/data/median_mdc.tsv'
# features = ['BALANCE_path_length', 'CHOREA-DOMINANT_path_length', 'CHOREA-NONDOMINANT_path_length', 'DAS-DOMINANT_SPIRAL_cv_vel', 'DAS-NONDOMINANT_SPIRAL_cv_vel', 'STT-DOMINANT_intertap_interval_mean', 'STT-NONDOMINANT_intertap_interval_mean', 'TMWT_step_freq_variance', 'UTT_step_freq_variance']

data = pd.read_csv(input_file, sep='\t', dtype={'subject_id': str},
                   usecols=['subject_id', 'weeks', 'feature', 'change_from_baseline', 
                            'arm', 'age12cap1'])

data = data[data['feature'].isin(features)] \
                           .dropna() \
                           .sort_values(by=['subject_id', 'feature', 'weeks'], ascending=True)

mdc = pd.read_csv(mdc_file, sep='\t', usecols=['test', 'feature_name', 'median_MDC'])

mdc['feature'] = mdc['test'] + '_' + mdc['feature_name']
merged = data.merge(mdc[['feature','median_MDC']], on=['feature'], how='left')

merged = merged[~(merged['weeks']==0)]
merged['event'] = merged['change_from_baseline'] <= merged['median_MDC']
merged['event'] = merged['event'].replace({True: 1, False: 0})

result = merged.groupby(['subject_id', 'feature', 'arm', 'age12cap1']) \
               .apply(time_to_event) \
               .reset_index()
    
result.to_csv(output_data, sep='\t', index=False)
