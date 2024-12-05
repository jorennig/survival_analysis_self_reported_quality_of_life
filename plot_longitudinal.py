import pandas as pd
from plotnine import ggplot, aes, geom_line, geom_errorbar, geom_point, xlab, ylab, facet_wrap, scale_color_manual

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]
feature = snakemake.wildcards.feature
treatment_arm = snakemake.wildcards.treatment_arm
flip_sign = snakemake.params.flip_sign[feature]
data = snakemake.wildcards.data
selection_criteria = snakemake.params.selection_criteria[data]

data = pd.read_csv(input_file, dtype = {'subject_id': str})
data = data[data['feature_name']==feature]

data = data[(data['arm']==treatment_arm) | (data['arm']=='PBO')]

data['select_data'] = data['change_from_baseline'] * flip_sign
data = data[data['select_data']>=selection_criteria]

summary = data.groupby(['feature_name', 'weeks', 'arm', 'age12cap1'])['change_from_baseline'] \
              .agg(['mean', 'median', 'std', 'sem', 'count', 'mad']) \
              .reset_index()
summary['y_lim_up'] = summary['mean'] + summary['sem']
summary['y_lim_low'] = summary['mean'] - summary['sem']

summary = summary.sort_values(by=['age12cap1'], ascending=False)
summary['age12cap1'] = pd.Categorical(summary['age12cap1'], ordered=True, 
                                      categories=summary['age12cap1'].unique())

width = 20
height = 5
plot = ggplot(summary, aes(x='weeks', y='mean', color='arm')) \
       + geom_line() \
       + geom_point() \
       + geom_errorbar(aes(ymin='y_lim_low', ymax='y_lim_up'), width=.2) \
       + xlab('Weeks') \
       + ylab('Change from Baseline') \
       + facet_wrap('age12cap1', ncol=summary['age12cap1'].nunique()) \
       + scale_color_manual(values=['royalblue', 'darkorange'])
plot

plot.save(output_file, dpi=300, width=width, height=height, verbose=False)
