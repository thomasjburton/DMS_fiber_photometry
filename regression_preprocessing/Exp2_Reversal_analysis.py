# %%
#Thomas Burton (t.burton@unsw.edu.au) Jan 2023
#Keep this file in .\scripts in your project folder e.g. .\Exp2_Reversal_BIDS\scripts
import logging
from pathlib import Path
import pandas as pd
import numpy as np
from behapy.utils import load_preprocessed_experiment
from behapy.events import build_design_matrix, regress, find_events
import statsmodels.api as sm
import seaborn as sns
sns.set_theme()


logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

 

# %%
BIDSROOT = Path('..')
pre = load_preprocessed_experiment(BIDSROOT)
dff_id = ['subject', 'session', 'task', 'run', 'label']
dff_recordings = pre.recordings.loc[:, dff_id].drop_duplicates()
#remove excluded animals
dff_recordings = dff_recordings.loc[~dff_recordings.subject.isin(['3', '5', '10', '12'])]

# %%
# z-score the dff
dff = pre.dff.copy()
dff['dff'] = dff.dff.groupby(dff_id, group_keys=False).apply(lambda df: (df - df.mean()) / df.std())


# %%
# Map IPSI and CONTRA
def _map_ipsi_contra(row):
    r = row.iloc[0]
    if r.label == 'RDMS':
        events = pre.events.loc[(r.subject, r.session, r.task, r.run, r.label)].replace({'rlp': 'ipsilp', 'llp': 'contralp'})
    elif r.label == 'LDMS':
        events = pre.events.loc[(r.subject, r.session, r.task, r.run, r.label)].replace({'rlp': 'contralp', 'llp': 'ipsilp'})
    else:
        raise ValueError(f'Unknown label {r.label}')
    return events.sort_index(level='onset')


events = dff_recordings.groupby(dff_id).apply(_map_ipsi_contra)
events

# %%
# Map events to individual recordings
def _get_nonevent(events, sub_events):
    nonevent = events.loc[:, ['duration']].merge(sub_events.loc[:, ['latency']], how='left', left_index=True, right_index=True, indicator=True)
    return nonevent.loc[nonevent._merge == 'left_only', ['duration', 'latency']]

REWmag = find_events(events, 'mag', ['pel', 'suc'])
Pelmag = find_events(events, 'mag', 'pel')
Sucmag = find_events(events, 'mag', 'suc')
NOREWmag = _get_nonevent(events.loc[events.event_id == 'mag', :], REWmag)
new_events = pd.concat([Pelmag, Sucmag, REWmag, NOREWmag],
                       keys=['Pelmag', 'Sucmag', 'REWmag', 'NOREWmag'],
                       names=['event_id'])
new_events = new_events.reset_index('event_id').loc[:, ['duration', 'event_id']]
events = pd.concat([events, new_events]).sort_index()
rew = events.loc[events.event_id.isin(['pel', 'suc'])].copy()
rew.event_id = 'rew'
events = pd.concat([events, rew]).sort_index()


# %%
plot_meta = {'Magazine': ['REWmag', 'NOREWmag'],
              'Reward': ['rew'],
              'Lever press': ['ipsilp', 'contralp']}

event_ids_of_interest = sum(plot_meta.values(), [])
events_of_interest = events.loc[events.event_id.isin(event_ids_of_interest), :]

# %%
def _build_design_matrix(row):
    r = row.iloc[0]
    return build_design_matrix(
        dff.loc[(r.subject, r.session, r.task, r.run, r.label), :],
        events_of_interest.loc[(r.subject, r.session, r.task, r.run, r.label), :],
        (-2.5, 2.5))


design_matrix = dff_recordings.groupby(dff_id).apply(_build_design_matrix).fillna(False).astype(bool)

# %%
idx = pd.IndexSlice
#Example of design matrix that extracts data from the three training tasks: FI15, RR5 and RR10 
design_filt = design_matrix.loc[idx[:, :, ['FI15', 'RR5', 'RR10'], :, :, :], :].sort_index()

#Example of design matrix that extracts the Reversal task data
#design_filt = design_matrix.loc[idx[:, :, ['Rev1'], :, :, :], :].sort_index()

def _regress(df):
    return regress(df, dff.loc[df.index, 'dff'], min_events=25)


# %%
#Several examples provided for grouping data

# Group data by sub only
#r1 = design_filt.loc[:, idx[sum(plot_meta.values(), []), :]].groupby(level=('subject'), group_keys=True).apply(_regress)

# Group data by sub and task
r1 = design_filt.loc[:, idx[sum(plot_meta.values(), []), :]].groupby(level=('subject', 'task'), group_keys=True).apply(_regress)

#Group data by sub, session and task
#r1 = design_filt.loc[:, idx[sum(plot_meta.values(), []), :]].groupby(level=('subject', 'session','task'), group_keys=True).apply(_regress)

# %%
#Toggle between the two s1 lines below depending on how data is structured in r1
s1 = r1.stack(0).stack()    
#s1 = r1.copy()
s1.name = 'beta'
s1 = s1.reset_index()
s1['event_type'] = s1.event.map({v: k for k, l in plot_meta.items() for v in l})
sns.relplot(data=s1, x='offset', y='beta', hue='event', row='event_type', col='task',
           kind='line', hue_order=sum(plot_meta.values(), []), aspect=2)


# %%

#Example of how to unstack and export data to a csv for further handling in MATLAB. Will need to adjust depending on s1 structure
r1.stack(0).to_csv('training_beta') 


