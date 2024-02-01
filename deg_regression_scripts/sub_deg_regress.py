# %%
from config import detrend,\
                   gen_session_mask,\
                   load_session_meta,\
                   load_fibre_session, load_events,\
                   gen_session_mask, find_events_adjacent,\
                   find_events_within,\
                   generate_all_event_regressors
import logging
from pathlib import Path
import pandas as pd
from sklearn.linear_model import LinearRegression

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()
%matplotlib inline

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
idx = pd.IndexSlice

# %%
PROJECTPATH = Path('..')
DSPATH = PROJECTPATH / 'derivatives/downsampled'

# %%
# Load sessions and cutoff info
sessions, cb_info = load_session_meta(PROJECTPATH/'etc/session-info.csv',
                                      PROJECTPATH/'etc/session-cutoffs.csv',
                                      PROJECTPATH/'etc/counterbalance-info.csv')

# Filter sessions for phases we're interested in
sessions_deg = sessions.loc[idx[:, :, ['DEG'], :, :]]

# %%
def regress_all_events(signal_df, events_df, mask, fs, window=(-1, 2)):
    df = generate_all_event_regressors(signal_df, events_df, fs, window)
    df.loc[mask] = False
    filt_df = df.loc[df.any(axis=1)]
    lr = LinearRegression().fit(filt_df.to_numpy(),
                                signal_df.loc[filt_df.index].to_numpy())
    return pd.Series(lr.coef_, index=df.columns)


# %%
def code_events(index, mag, pel, suc, rlp, llp, ipsi_rew, contra_rew,
                deg_rew, nondeg_rew,
                ipsi_lp, contra_lp, deg_lp, nondeg_lp):
    df = pd.DataFrame(index=index)
    rew = pel | suc
    lp = rlp | llp

    df['ipsi_deg_lp'] = ipsi_lp & deg_lp
    df['contra_deg_lp'] = contra_lp & deg_lp
    df['ipsi_nondeg_lp'] = ipsi_lp & nondeg_lp
    df['contra_nondeg_lp'] = contra_lp & nondeg_lp

    df['DEGREWmag'] = find_events_adjacent(deg_rew, mag)
    df['NONDEGREWmag'] = find_events_adjacent(nondeg_rew, mag)

    df['DEGLPmag'] = find_events_adjacent(deg_lp, mag)
    df['NONDEGLPmag'] = find_events_adjacent(nondeg_lp, mag)
    df.loc[df.DEGLPmag & df.NONDEGLPmag,
           ['DEGLPmag', 'NONDEGLPmag']] = False
    df['EARNEDrew'] = find_events_within(lp, rew, tmax=1., allow_simultaneous=True)
    df['FREErew'] = rew & ~df.EARNEDrew

    df['REWmag'] = find_events_adjacent(rew, mag)
    df['EMPTYmag'] = mag & ~df.REWmag
    df['EARNEDREWmag'] = find_events_adjacent(df.EARNEDrew, mag)
    df['FREEREWmag'] = find_events_adjacent(df.FREErew, mag)

    df['EARNEDDEGREWmag'] = df.EARNEDREWmag & df.DEGREWmag
    df['EARNEDNONDEGREWmag'] = df.EARNEDREWmag & df.NONDEGREWmag
    df['IPSILPmag'] = find_events_adjacent(ipsi_lp, mag)
    df['CONTRALPmag'] = find_events_adjacent(contra_lp, mag)

    return df


# %%
def load_all_sessions(sessions, cb_info, path, window=(-1., 2.)):
    dffs = []
    events = []
    ev_regs = []
    fs = None
    for info in sessions.itertuples():
        ses_info = (info.animal_id, info.session_id, info.task, info.task_id)
        print(ses_info)
        mask_info = (info.t_start, info.t_end,
                     info.lower_bound, info.upper_bound)
        df, fs = load_fibre_session(path, *ses_info)
        mask = gen_session_mask(df.loc['405', :], *mask_info)
        mask.name = 'mask'
        dff = pd.Series(df.loc['465', :] / df.loc['465', :].mean(),
                        index=df.loc['465', :].index, dtype=float, name='dff')
        dff.loc[:] = detrend(dff, fs)
        dff = (dff / dff.std()).to_frame()
        dff['mask'] = mask
        dffs.append(dff)
        ev = load_events(path, *ses_info)
        sub_cb = cb_info.loc[cb_info.animal_id == info.animal_id]
        ev_coded = code_events(ev.index, ev.mag, ev.pel, ev.suc, ev.rlp,
                               ev.llp,
                               ev[sub_cb.ipsi_rew.iloc[0]],
                               ev[sub_cb.contra_rew.iloc[0]],
                               ev[sub_cb.deg_rew.iloc[0]],
                               ev[sub_cb.nondeg_rew.iloc[0]],
                               ev[sub_cb.ipsi_lp.iloc[0]],
                               ev[sub_cb.contra_lp.iloc[0]],
                               ev[sub_cb.deg_lp.iloc[0]],
                               ev[sub_cb.nondeg_lp.iloc[0]])
        events.append(ev_coded)
        ev_reg = generate_all_event_regressors(dff.dff, ev_coded, fs, window=window)
        ev_reg.loc[dff['mask']] = False
        ev_reg = ev_reg.loc[ev_reg.any(axis=1)]
        ev_regs.append(ev_reg)
    dff_df = pd.concat(dffs, axis=0,
                       keys=zip(sessions.animal_id, sessions.session_id,
                                sessions.task, sessions.task_id,
                                sessions.tank),
                       names=('animal_id', 'session_id', 'task', 'task_id',
                              'tank_id'))
    events_df = pd.concat(events, axis=0,
                          keys=zip(sessions.animal_id, sessions.session_id,
                                   sessions.task, sessions.task_id,
                                   sessions.tank),
                          names=('animal_id', 'session_id', 'task', 'task_id',
                                 'tank_id'))
    regs_df = pd.concat(ev_regs, axis=0,
                        keys=zip(sessions.animal_id, sessions.session_id,
                                 sessions.task, sessions.task_id,
                                 sessions.tank),
                        names=('animal_id', 'session_id', 'task', 'task_id',
                               'tank_id'))
    return dff_df, events_df, regs_df, fs


# %%
def regress_mass(df, dff_df, min_events=100):
    _df = df.loc[:, df.sum() > min_events]
    if _df.empty:
        return pd.Series(dtype=float, index=_df.columns)
    lr = LinearRegression().fit(_df.to_numpy(),
                                dff_df.dff.loc[df.index].to_numpy())
    return pd.Series(lr.coef_, index=_df.columns)


# %%
# Load contingency degradation sessions
dff_df, events_df, regressors_df, fs = \
    load_all_sessions(sessions_deg.reset_index(), cb_info, DSPATH, window=(-2.5, 2.5))

# %%
fig_path = PROJECTPATH / 'derivatives/figures'

# %%
# Let's look at press types, mag checks by press type and mag checks by
# reward, between subjects across 'session_group' (groups of tanks).
plot_meta = {'Lateral magazine': ['IPSILPmag', 'CONTRALPmag'],
             'Magazine': ['EARNEDDEGREWmag', 'EARNEDNONDEGREWmag',
                          'FREEREWmag', 'EMPTYmag'],
             'Lever press': ['ipsi_deg_lp', 'contra_deg_lp',
                             'ipsi_nondeg_lp', 'contra_nondeg_lp']}
r1 = regressors_df.loc[idx[:, :, :, :, :, :], idx[sum(plot_meta.values(), []), :]].groupby(level=('animal_id'), group_keys=True).apply(regress_mass, dff_df=dff_df, min_events=25)
# s1 = r1.stack(0).stack()
s1 = r1.copy()
s1.name = 'Regression coeffs'
s1 = s1.reset_index()
s1['event_type'] = s1.event.map({v: k for k, l in plot_meta.items() for v in l})
fig = sns.relplot(data=s1, x='offset', y='Regression coeffs', hue='event',
                  row='event_type', kind='line',
                  hue_order=sum(plot_meta.values(), []))

fig.savefig(fig_path / 'deg_massregress_lplatbydeg-maglat-magtype_alltanks.png', dpi=300, format='png')
r1out = pd.concat([r1.groupby(['event', 'offset']).mean(),
                   r1.groupby(['event', 'offset']).sem()],
                  axis=1, keys=['mean', 'sem'], names=['stats'])
r1out = r1out.stack().unstack(1)
r1out.to_csv(fig_path / 'deg_massregress_lplatbydeg-maglat-magtype_alltanks_wide.csv')
r1.unstack().to_csv(fig_path / 'deg_massregress_lplatbydeg-maglat-magtype_alltanks_subjects.csv')
