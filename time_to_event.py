import pandas as pd

def time_to_event(df):    
    df['shift'] = df['change_from_baseline'].shift(1)    
    df['event'] = (df['shift']>=df['change_from_baseline']) & \
                  (df['change_from_baseline']>=df['change'])    
    event_detected = df['event'].any()
    weeks_to_event = df.loc[df['event'],'weeks'].min()
    return pd.Series([weeks_to_event, event_detected], index=['weeks_to_event', 'event_detected'])

input_file = snakemake.input.data[0]
output_data = snakemake.output[0]
feature = snakemake.wildcards.feature
flip_sign = snakemake.params.flip_sign[feature]
change_value = snakemake.params.change_value[feature]

data = pd.read_csv(input_file, dtype={'subject_id': str})
data = data[data['feature_name']==feature]

data['change'] = change_value
data['change_from_baseline'] = data['change_from_baseline'] * flip_sign

result = data.groupby(['subject_id', 'feature_name', 'arm', 'age12cap1']) \
             .apply(time_to_event) \
             .reset_index()
result['weeks_to_event'] = result['weeks_to_event'] - 1
result['weeks_to_event'] = result['weeks_to_event'].fillna(data['weeks'].max())
result['all_patients'] = 'all_patients'

result.to_csv(output_data, index=False)
