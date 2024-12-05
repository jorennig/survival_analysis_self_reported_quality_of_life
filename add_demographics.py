import pandas as pd

input_file = snakemake.input.data
input_demo = snakemake.input.demo
output_data = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})

demo = pd.read_csv(input_demo, dtype={'subject_id': str}) \
         .filter(['subject_id', 'arm', 'age12cap1'])

data = data.sort_values(by=['subject_id', 'feature_name', 'two_week_period'])

data = data.merge(demo, on='subject_id')

data.to_csv(output_data, index=False)
