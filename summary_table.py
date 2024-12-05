import pandas as pd

input_file = r'\\exports.hps.kau.science.roche.com\pred\rpmda\BN40423-v2\jira\RPMDA-8747-eq5d5l-survival-analysis\data\eq5d5l_change_from_baseline.csv'

data = pd.read_csv(input_file, dtype={'subject_id': str})

data = data[data['weeks']==1]

data = data[data['feature_name'].isin(['VAS_score', 'index_score'])]

summary = data.groupby(['feature_name', 'arm', 'age12cap1'])['baseline'].mean().reset_index(name='mean_baseline')
summary = pd.pivot_table(summary, values='mean_baseline', index=['feature_name', 'arm'], columns=['age12cap1']).reset_index()
