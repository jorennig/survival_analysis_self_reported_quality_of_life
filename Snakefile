import pandas as pd
import os

configfile: "config.yaml"

rule add_demographics:
    input: script='scripts/add_demographics.py',
           data=config['data_file'],
           demo=config['demographics']
    output: 'data/eq5d5l_demographics.csv'
    script: 'scripts/add_demographics.py'

rule format:
    input: script='scripts/format.py',
           data=rules.add_demographics.output
    output: 'data/eq5d5l_format.csv'
    script: 'scripts/format.py'

rule change_from_baseline:
    input: script='scripts/change_from_baseline.py',
           data=rules.format.output
    output: 'data/eq5d5l_change_from_baseline.csv'
    params: max_week=config['max_week']
    script: 'scripts/change_from_baseline.py'

rule time_to_event:
    input: script='scripts/time_to_event.py',
           data=rules.change_from_baseline.output
    output: 'data/{feature}_time_to_event.csv'
    params: flip_sign=config['flip_sign'], 
            change_value=config['change_value']
    script: 'scripts/time_to_event.py'

rule concatenate:
    input: expand('data/{feature}_time_to_event.csv', feature=config['features'])
    output: 'data/time_to_event_concatenated.csv'
    run: pd.concat([pd.read_csv(f, dtype={'subject_id': str}) for f in input]).to_csv(output[0], index=False)

rule summary_events:
    input: script='scripts/summary_events.py',
           data=rules.concatenate.output
    output: 'results/{analysis}_summary_events.csv'
    params: grouping=config['grouping']
    script: 'scripts/summary_events.py'

rule log_rank_test:
    input: script='scripts/log_rank_test.py',
           data=rules.time_to_event.output
    output: 'results/{treatment_arm}/{analysis}/{feature}_log_rank_test.csv'
    params: grouping=config['grouping'],
            control_arm=config['control_arm']
    script: 'scripts/log_rank_test.py'

rule kaplan_meier_plot_all_patients:
    input: script='scripts/kaplan_meier_plot.py',
           data=rules.time_to_event.output,
           results='results/{treatment_arm}/all_patients/{feature}_log_rank_test.csv'
    output: 'results/{treatment_arm}/all_patients/{feature}_kaplan_meier_plot.png'
    params: control_arm=config['control_arm'],
            group_var='all_patients'
    script: 'scripts/kaplan_meier_plot.py'

rule kaplan_meier_plot_disease_group:
    input: script='scripts/kaplan_meier_plot.py',
           data=rules.time_to_event.output,
           results='results/{treatment_arm}/disease_group/{feature}_log_rank_test.csv'
    output: 'results/{treatment_arm}/disease_group/{feature}_kaplan_meier_plot.png'
    params: control_arm=config['control_arm'],
            group_var='age12cap1'
    script: 'scripts/kaplan_meier_plot.py'

rule plot_longitudinal:
    input: script='scripts/plot_longitudinal.py',
           data=rules.change_from_baseline.output,
    output: 'results/longitudinal_{treatment_arm}_{data}/{feature}_longitudinal.png'
    params: flip_sign=config['flip_sign'],
            selection_criteria=config['selection_criteria']
    script: 'scripts/plot_longitudinal.py'


rule all:
    input: expand('results/{analysis}_summary_events.csv', analysis=config['analyses']),
           expand('results/{treatment_arm}/all_patients/{feature}_kaplan_meier_plot.png', feature=config['features'], treatment_arm=config['treatment_arms']),
           expand('results/{treatment_arm}/disease_group/{feature}_kaplan_meier_plot.png', feature=config['features'], treatment_arm=config['treatment_arms']),
           expand('results/longitudinal_{treatment_arm}_{data}/{feature}_longitudinal.png', feature=config['features'], data=config['data_selection'], treatment_arm=config['treatment_arms'])
