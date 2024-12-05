import pandas as pd

input_file = snakemake.input.data[0]
output_data = snakemake.output[0]
max_week = int(snakemake.params.max_week)

data = pd.read_csv(input_file, dtype={'subject_id': str})

baseline = data[data['weeks']==0] \
               .filter(['subject_id', 'feature_name', 'numeric_value']) \
               .rename(columns={'numeric_value': 'baseline'}) \
               .drop_duplicates(['subject_id', 'feature_name'], keep='first')

data = data.merge(baseline, on=['subject_id', 'feature_name'])

data['change_from_baseline'] = data['numeric_value'] - data['baseline']

data = data[(data['weeks']>0) & (data['weeks']<max_week)]

data.to_csv(output_data, index=False)
