import pandas as pd
import numpy as np
from lifelines.statistics import logrank_test

def log_rank(chunk, c, t):
    pcb = chunk[chunk['arm'] == c]
    treatment = chunk[chunk['arm'] == t]
    
    pcb_weeks = np.array(pcb['weeks_to_event'])
    pcb_event = np.array(pcb['event_detected'])    
    treatment_weeks = np.array(treatment['weeks_to_event'])
    treatment_event = np.array(treatment['event_detected'])
    
    lr_test = logrank_test(pcb_weeks, treatment_weeks, pcb_event, treatment_event)
    return pd.Series([t, lr_test.p_value], index=['treatment_arm', 'p_value'])

input_file = snakemake.input.data[0]
output_data = snakemake.output[0]
analysis = snakemake.wildcards.analysis
grouping = snakemake.params.grouping[analysis]
treatment_arm = snakemake.wildcards.treatment_arm
control_arm = snakemake.params.control_arm

dtype={'subject_id': str, 'weeks_to_event': int, 'event_detected': int}
data = pd.read_csv(input_file, dtype=dtype)

results = data.groupby(grouping).apply(log_rank, control_arm, treatment_arm).reset_index()
results['all_patients'] = 'all_patients'

results.to_csv(output_data, index=False)
