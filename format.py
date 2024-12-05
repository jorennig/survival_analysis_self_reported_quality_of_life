import pandas as pd

input_file = snakemake.input.data[0]
output_data = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})

age12cap1_dict = {1: 'low_age_low_cap', 2: 'low_age_high_cap',
                  3: 'high_age_low_cap', 4: 'high_age_high_cap'}
data['age12cap1'] = data['age12cap1'].replace(age12cap1_dict)

replace = '|'.join(['/', ' ', '-'])
data['feature_name'] = data['feature_name'].str.replace(replace, '_')

data['weeks'] = data['day_of_study'] // 7

data.to_csv(output_data, index=False)
