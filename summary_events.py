import pandas as pd

input_file = snakemake.input.data[0]
output_data = snakemake.output[0]
analysis = snakemake.wildcards.analysis
grouping = snakemake.params.grouping[analysis]

data = pd.read_csv(input_file, dtype={'subject_id': str})

n = data.groupby(grouping)['subject_id'] \
        .nunique() \
        .reset_index(name='n')

summary = data.groupby(grouping)['event_detected'] \
              .sum() \
              .reset_index()

summary = summary.merge(n, on=grouping)

summary.to_csv(output_data, index=False)
