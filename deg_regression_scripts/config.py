import json
import numpy as np
import pandas as pd
import scipy.signal as sig


def load_session_info(filename):
    dtypes = {
        'Path': str,
        'AnimalA': pd.Int32Dtype(),
        'AnimalB': pd.Int32Dtype(),
        'Tank': pd.Int32Dtype(),
        'Task': str,
        'SessionA': pd.Int32Dtype(),
        'SessionB': pd.Int32Dtype()
    }
    info = pd.read_csv(filename, dtype=dtypes)
    column_map = {'Path': 'path',
                  'AnimalA': 'animal_id',
                  'AnimalB': 'animal_id',
                  'Tank': 'tank',
                  'Task': 'task',
                  'SessionA': 'task_id',
                  'SessionB': 'task_id'}
    longform = pd.concat([
        info[['Path', 'AnimalA', 'Tank', 'Task', 'SessionA']].rename(column_map, axis=1),
        info[['Path', 'AnimalB', 'Tank', 'Task', 'SessionB']].rename(column_map, axis=1)],
        keys=['A', 'B'], names=['port', 'csv_row'])
    longform['session_id'] = longform.groupby('animal_id').cumcount()
    return longform


def load_session_cutoffs(filename):
    dtypes = {
        'Path': str,
        'PrtA_lowerBound': pd.Int32Dtype(),
        'PrtB_lowerBound': pd.Int32Dtype(),
        'PrtA_upperBound': pd.Int32Dtype(),
        'PrtB_upperBound': pd.Int32Dtype(),
        'PrtA_timeCutoff': pd.Int32Dtype(),
        'PrtB_timeCutoff': pd.Int32Dtype()
    }
    column_map = {
        'Path': 'path',
        'PrtA_lowerBound': 'lower_bound',
        'PrtB_lowerBound': 'lower_bound',
        'PrtA_upperBound': 'upper_bound',
        'PrtB_upperBound': 'upper_bound',
        'PrtA_timeCutoff': 't_end',
        'PrtB_timeCutoff': 't_end'
    }
    cutoffs = pd.read_csv(filename, dtype=dtypes)
    longform = pd.concat([
        cutoffs[['Path', 'PrtA_lowerBound', 'PrtA_upperBound', 'PrtA_timeCutoff']].rename(column_map, axis=1),
        cutoffs[['Path', 'PrtB_lowerBound', 'PrtB_upperBound', 'PrtB_timeCutoff']].rename(column_map, axis=1)],
        keys=['A', 'B'], names=['port', 'csv_row'])
    longform['t_start'] = 8
    return longform


def load_counterbalance_info(filename):
    dtypes = {
        'animal_id': pd.Int32Dtype(),
        'port': str,
        'run': pd.Int32Dtype(),
        'deg_outcome': str,
        'deg_lever': str,
        'pel_lever': str,
        'canula_side': str
    }
    df = pd.read_csv(filename, dtype=dtypes)
    df['ipsi_rew'] = np.where(df.pel_lever == df.cannula_side, 'pel', 'suc')
    df['contra_rew'] = np.where(df.pel_lever == df.cannula_side, 'suc', 'pel')
    df['ipsi_lp'] = np.where(df.cannula_side == 'L', 'llp', 'rlp')
    df['contra_lp'] = np.where(df.cannula_side == 'L', 'rlp', 'llp')
    df['deg_rew'] = np.where(df.deg_outcome == 'P', 'pel', 'suc')
    df['nondeg_rew'] = np.where(df.deg_outcome == 'P', 'suc', 'pel')
    df['deg_lp'] = np.where(df.deg_lever == 'L', 'llp', 'rlp')
    df['nondeg_lp'] = np.where(df.deg_lever == 'L', 'rlp', 'llp')
    return df


def detrend(x, fs):
    try:
        b = detrend.filter_b
    except AttributeError:
        b = sig.firwin2(501, freq=[0, 0.05, 0.1, fs/2],
                        gain=[0., 0.001, 1., 1.], fs=fs)
        detrend.filter_b = b
    return sig.filtfilt(b, 1, x)


def load_session_meta(session_fn, cutoffs_fn, counterbalance_fn):
    info = load_session_info(session_fn)
    cutoffs = load_session_cutoffs(cutoffs_fn)
    cb_info = load_counterbalance_info(counterbalance_fn)
    cutoff_columns = ['lower_bound', 'upper_bound', 't_start', 't_end']
    info = info.join(cutoffs[cutoff_columns])
    session_fields = ['animal_id', 'tank', 'task', 'task_id', 'session_id',
                    'lower_bound', 'upper_bound', 't_start', 't_end']
    sessions = info.loc[info.animal_id != 0, session_fields]
    index_columns = ['animal_id', 'session_id', 'task', 'task_id', 'tank']
    return sessions.set_index(index_columns), cb_info


def gen_session_mask(df, t_start, t_end, lower_bound, upper_bound):
    mask = df.index.to_frame().time.between(t_start, t_end)
    mask &= df.between(lower_bound, upper_bound)
    return ~mask


def load_fibre_session(base, subject, session, task, task_id):
    ses_path = base / 'sub-{sub:02d}/ses-{ses:02d}'.format(sub=subject, ses=session)
    fp_template = 'fp/sub-{sub:02d}_ses-{ses:02d}_task-{task}-{task_id:02d}_L{ch}.{ext}'
    series = []
    channels = ['405', '465']
    fs = -1.
    for ch in channels:
        meta_fn = fp_template.format(sub=subject, ses=session, task=task,
                                     task_id=task_id, ch=ch, ext='json')
        data_fn = fp_template.format(sub=subject, ses=session, task=task,
                                     task_id=task_id, ch=ch, ext='npy')
        with open(ses_path/meta_fn) as file:
            meta = json.load(file)
        data = np.load(ses_path/data_fn)
        if (fs > 0) and (meta['fs'] != fs):
            logging.warning('Channels have different frequencies')
        fs = meta['fs']
        ts = np.arange(data.shape[0]) / fs
        series.append(pd.Series(data, index=ts))
    return (pd.concat(series, axis=0, keys=channels, names=['channel', 'time']),
            fs)


def extract_events(event_ints):
    event_bits = [list('{:08b}'.format(x)) for x in event_ints]
    return np.array(event_bits)


def find_events_adjacent(event1_mask, event2_mask, direction='forward'):
    df1 = pd.DataFrame(index=np.where(event1_mask)[0])
    df2 = pd.DataFrame(np.where(event2_mask)[0], index=np.where(event2_mask)[0],
                       columns=['indices'])
    indices = pd.merge_asof(df1, df2,
                            left_index=True, right_index=True,
                            direction=direction).indices
    series = pd.Series([False]*len(event1_mask))
    series.loc[indices[~indices.isna()]] = True
    return series.to_numpy()


def find_events_within(event1_mask: pd.Series, event2_mask: pd.Series,
                       direction='forward', tmax=1.,
                       allow_simultaneous=False):
    # Inputs need to be Series with timestamps
    event1_mask.name = 'src'
    event2_mask.name = 'dest'
    df1 = pd.DataFrame(index=event1_mask.index[event1_mask])
    df1.index.name = 'onset'
    # Subtract or add one microsecond to ensure the same event isn't
    # detected as following or preceding itself (if the same event train
    # is passed in)
    if allow_simultaneous:
        offset = 0
    else:
        offset = -1e-6 if direction == 'forward' else 1e-6
    df2 = pd.DataFrame(event2_mask.index[event2_mask].to_numpy(),
                       index=(event2_mask.index[event2_mask] + offset),
                       columns=['indices'])
    merged = pd.merge_asof(df1, df2,
                           left_index=True, right_index=True,
                           direction=direction).reset_index()
    merged = merged[(merged.indices - merged.onset).abs() < tmax]
    series = pd.Series([False]*len(event2_mask), index=event2_mask.index)
    series.loc[merged.indices[~merged.indices.isna()]] = True
    return series.to_numpy()


def find_events_with_adjacent(event1_mask, event2_mask, direction='forward'):
    df1 = pd.DataFrame(index=np.where(event1_mask)[0])
    df2 = pd.DataFrame(np.where(event2_mask)[0], index=np.where(event2_mask)[0],
                       columns=['indices'])
    index = pd.merge_asof(df1, df2,
                          left_index=True, right_index=True,
                          direction=direction).index
    series = pd.Series([False]*len(event1_mask))
    series.loc[index[~index.isna()]] = True
    return series.to_numpy()


def load_events(base, subject, session, task, task_id,
                event_names=['ue0', 'ue1', 'ue2', 'pel', 'suc', 'rlp', 'mag', 'llp']):
    ses_path = base / 'sub-{sub:02d}/ses-{ses:02d}'.format(sub=subject, ses=session)
    events_template = 'sub-{sub:02d}_ses-{ses:02d}_task-{task}-{task_id:02d}_events.tsv'
    events_fn = ses_path/events_template.format(sub=subject, ses=session,
                                                task=task, task_id=task_id)
    df = pd.read_csv(events_fn, delimiter='\t', index_col=0, dtype=float)
    df['value'] = df['value'].astype(int)
    df = df.reset_index(drop=True)
    events = extract_events(df['value']).astype(int)
    event_onsets = np.concatenate([events[0:1, :], np.diff(events, axis=0)]) == 1
    onsets_df = pd.DataFrame(event_onsets,
                             columns=event_names,
                             index=df.onset)
    return onsets_df[event_names[:]]
    # return onsets_df[event_names[3:]]


def find_nearest(origin, fit):
    df0 = pd.DataFrame(np.array(origin), index=origin, columns=['origin'])
    df1 = pd.DataFrame(np.array(fit), index=fit, columns=['fit'])
    first = pd.merge_asof(df1, df0,
                          left_index=True, right_index=True,
                          direction='nearest')
    second = pd.DataFrame(index=origin, dtype=bool)
    second['nearest'] = False
    second.loc[first['origin'], 'nearest'] = True
    return second['nearest'].to_numpy()


def generate_event_regressors(event_mask, fs, window=(-1, 2)):
    # Make sure event_mask is an array
    event_mask = np.array(event_mask)
    window_indices = np.round(np.array(window) * fs).astype(int)
    regressor_offsets = np.arange(*window_indices)
    regressors = np.zeros(
        np.concatenate([event_mask.shape, regressor_offsets.shape]))
    for i, offset in enumerate(regressor_offsets):
        if offset < 0:
            regressors[:offset, i] = event_mask[-offset:]
        elif offset > 0:
            regressors[offset:, i] = event_mask[:-offset]
        else:
            regressors[:, i] = event_mask
    return regressors, regressor_offsets


def generate_all_event_regressors(signal_df, events_df, fs, window=(-1, 2)):
    # Make this not accept dataframes for flexibility
    regressor_dfs = []
    for event in events_df:
        nearest = find_nearest(signal_df.index,
                               events_df.index[events_df[event]])
        regressors, regressor_offsets = \
            generate_event_regressors(nearest, fs, window=window)
        column_index = pd.MultiIndex.from_product([[event], regressor_offsets/fs],
                                                  names=('event', 'offset'))
        _df = pd.DataFrame(regressors, dtype=bool, index=signal_df.index,
                           columns=column_index)
        regressor_dfs.append(_df)
    return pd.concat(regressor_dfs, axis=1)
