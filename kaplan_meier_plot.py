import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

input_file = snakemake.input.data[0]
result_file = snakemake.input.results
output_data = snakemake.output[0]
feature = snakemake.wildcards.feature
treatment_arm = snakemake.wildcards.treatment_arm
control_arm = snakemake.params.control_arm
group_var = snakemake.params.group_var

dtype={'subject_id': str, 'weeks_to_event': int, 'age12cap1': str, 'event_detected': int}
data = pd.read_csv(input_file, dtype=dtype)

results = pd.read_csv(result_file, dtype=dtype)

data = data.merge(results, on=group_var)

data = data.sort_values(by=group_var, ascending=False)

n = data[group_var].nunique()
plt.figure(figsize=(5*n, 5))

i=0
for group, chunk in data.groupby(group_var, sort=False):
    
    i=i+1
    pcb = chunk[chunk['arm'] == control_arm]
    treatment = chunk[chunk['arm'] == treatment_arm]
    
    pcb_weeks = np.array(pcb['weeks_to_event'])
    pcb_event = np.array(pcb['event_detected'])
    
    treatment_weeks = np.array(treatment['weeks_to_event'])
    treatment_event = np.array(treatment['event_detected'])
    
    ax = plt.subplot(1, n, i)
    kmf_pbo = KaplanMeierFitter().fit(pcb_weeks, event_observed=pcb_event, label=control_arm)
    ax = kmf_pbo.plot_survival_function(ax=ax)
    kmf_t = KaplanMeierFitter().fit(treatment_weeks, event_observed=treatment_event, label=treatment_arm)
    ax = kmf_t.plot_survival_function(ax=ax)
    ax.set_xlabel('weeks')
    plt.ylim(0.3, 1)
    p_value = chunk['p_value'].iloc[0]
    group = chunk[group_var].iloc[0]
    ax.set_title('%s\nlogrank p-value = %.4f' % (group, p_value))
    
plt.tight_layout()
plt.savefig(output_data, dpi=300, bbox_inches="tight")
plt.close()
