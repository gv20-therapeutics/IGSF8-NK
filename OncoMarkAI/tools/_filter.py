import pandas as pd


def sample_size_filter(df, groupby, cutoff):
    """
    Filter out patients with < cutoff samples

    Parameters
    ----------
    df : pd.DataFrame
        Data in "record" format. Each row represents a record. Indices represent groups.
    groupby : str, or list
        The key used to separate patients, must be present in the columns of the DataFrame
    cutoff : int
        The minimal number of samples
    Returns
    -------
    df : pd.DataFrame
        The filtered dataframe
    """

    tmp = df.groupby(groupby).size()
    tmp = tmp[tmp >= cutoff]
    output = df[df.set_index(groupby).index.isin(tmp.index)]

    return output
