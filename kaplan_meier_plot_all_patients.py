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

# input_file = r'\\exports.hps.kau.science.roche.com\pred\rpmda\BN40423-v2\jira\RPMDA-8747-eq5d5l-survival-analysis\data\Usual_activities_time_to_event.csv'
# result_file = r'\\exports.hps.kau.science.roche.com\pred\rpmda\BN40423-v2\jira\RPMDA-8747-eq5d5l-survival-analysis\results\Q16W\all_patients\Usual_activities_log_rank_test.csv'
# feature = 'Usual_activities'
# treatment_arm = 'Q16W'
# control_arm = 'PBO'

dtype={'subject_id': str, 'weeks_to_event': int, 'age12cap1': str, 'event_detected': int}
data = pd.read_csv(input_file, dtype=dtype)

results = pd.read_csv(result_file, dtype=dtype)

data = data.merge(results, on=['feature_name'])

pcb = data[data['arm'] == control_arm]
treatment = data[data['arm'] == treatment_arm]

pcb_weeks = np.array(pcb['weeks_to_event'])
pcb_event = np.array(pcb['event_detected'])

treatment_weeks = np.array(treatment['weeks_to_event'])
treatment_event = np.array(treatment['event_detected'])

plt.figure(figsize=(5, 5))

ax = plt.subplot(1, 1, 1)
kmf_pbo = KaplanMeierFitter().fit(pcb_weeks, event_observed=pcb_event, label=control_arm)
ax = kmf_pbo.plot_survival_function(ax=ax)
kmf_t = KaplanMeierFitter().fit(treatment_weeks, event_observed=treatment_event, label=treatment_arm)
ax = kmf_t.plot_survival_function(ax=ax)
ax.set_xlabel('weeks')
plt.ylim(0.3, 1)
p_value = data['p_value'].iloc[0]
ax.set_title('logrank p-value = %.4f' % (p_value))
    
plt.tight_layout()
plt.savefig(output_data, dpi=300, bbox_inches="tight")
plt.close()
