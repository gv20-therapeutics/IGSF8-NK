import numpy as np
import pandas as pd
from scipy.stats import (
    ttest_ind,
    mannwhitneyu,
    kruskal,
    f_oneway,
    spearmanr,
    pearsonr,
)
from statsmodels.stats.multitest import multipletests


def t_test(
    data: pd.DataFrame,
    groupby: str,
    group1: str,
    group2: str,
    value: str,
    min_sample_size1: int = 3,
    min_sample_size2: int = 3,
):
    """
    Compute the T-test for the means of two independent samples of scores.

    Parameters
    ----------
    data : pd.DataFrame
        Data in "record" format. Each row represents a record. Indices represent groups.
    groupby : str
        The name of the column that contains the group information
    group1 : str
        The name of the first group
    group2 : str
        The name of the second group
    min_sample_size1 : int
        The minimal sample size of the first group
    min_sample_size2 : int
        The minimal sample size of the second group
    value : str
        The name of the column that contains the value information

    Returns
    -------
    statistic
        The T statistic of the test
    pval
        The pval of the test
    mean
        The mean of each group
    std
        The standard deviation of each group
    median
        The median of each group
    sample_size
        The sample size of each group
    """

    # get the data
    data1 = data.query(f"`{groupby}` == '{group1}'")[value].dropna()
    data2 = data.query(f"`{groupby}` == '{group2}'")[value].dropna()

    # get the sample size
    sample_size = [len(data1), len(data2)]

    # calculate the mean
    mean = [np.mean(data1), np.mean(data2)]
    std = [np.std(data1), np.std(data2)]
    median = [np.median(data1), np.median(data2)]

    # perform the test
    if len(data1) >= min_sample_size1 and len(data2) >= min_sample_size2:
        statistic, pval = ttest_ind(data1, data2)
    else:
        statistic, pval = np.nan, np.nan

    # return the result
    return statistic, pval, mean, std, median, sample_size


def u_test(
    data: pd.DataFrame,
    groupby: str,
    group1: str,
    group2: str,
    value: str,
    min_sample_size1: int = 3,
    min_sample_size2: int = 3,
):
    """
    Compute the Mann-Whitney U test for two independent samples.

    Parameters
    ----------
    data : pd.DataFrame
        Data in "record" format. Each row represents a record. Indices represent groups.
    groupby : str
        The name of the column that contains the group information
    group1 : str
        The name of the first group
    group2 : str
        The name of the second group
    min_sample_size1 : int
        The minimal sample size of the first group
    min_sample_size2 : int
        The minimal sample size of the second group
    value : str
        The name of the column that contains the value information

    Returns
    -------
    statistic
        The U statistic of the test
    pval
        The pval of the test
    mean
        The mean of each group
    std
        The standard deviation of each group
    median
        The median of each group
    sample_size
        The sample size of each group
    """

    # get the data
    data1 = data.query(f"`{groupby}` == '{group1}'")[value].dropna()
    data2 = data.query(f"`{groupby}` == '{group2}'")[value].dropna()

    # get the sample size
    sample_size = [len(data1), len(data2)]

    # calculate the mean
    mean = [np.mean(data1), np.mean(data2)]
    std = [np.std(data1), np.std(data2)]
    median = [np.median(data1), np.median(data2)]

    # perform the test
    if len(data1) >= min_sample_size1 and len(data2) >= min_sample_size2:
        statistic, pval = mannwhitneyu(data1, data2)
    else:
        statistic, pval = np.nan, np.nan

    # return the result
    return statistic, pval, mean, std, median, sample_size


def pearson(
    data: pd.DataFrame,
    groupby: str,
    group1: str,
    group2: str,
    value: str,
    min_sample_size1: int = 3,
    min_sample_size2: int = 3,
):
    """
    Compute the pearson correlation for two independent samples.

    Parameters
    ----------
    data : pd.DataFrame
        Data in "record" format. Each row represents a record. Indices represent groups.
    groupby : str
        The name of the column that contains the group information
    group1 : str
        The name of the first group
    group2 : str
        The name of the second group
    min_sample_size1 : int
        The minimal sample size of the first group
    min_sample_size2 : int
        The minimal sample size of the second group
    value : str
        The name of the column that contains the value information

    Returns
    -------
    statistic
        The pearson correlation of the test
    pval
        The pval of the test
    mean
        The mean of each group
    std
        The standard deviation of each group
    median
        The median of each group
    sample_size
        The sample size of each group
    """

    # get the data
    data1 = data.query(f"`{groupby}` == '{group1}'")[value].dropna()
    data2 = data.query(f"`{groupby}` == '{group2}'")[value].dropna()

    # get the sample size
    sample_size = [len(data1), len(data2)]

    # calculate the mean
    mean = [np.mean(data1), np.mean(data2)]
    std = [np.std(data1), np.std(data2)]
    median = [np.median(data1), np.median(data2)]

    # perform the test
    if len(data1) >= min_sample_size1 and len(data2) >= min_sample_size2:
        statistic, pval = pearsonr(data1, data2)
    else:
        statistic, pval = np.nan, np.nan

    # return the result
    return statistic, pval, mean, std, median, sample_size


def spearman(
    data: pd.DataFrame,
    groupby: str,
    group1: str,
    group2: str,
    value: str,
    min_sample_size1: int = 3,
    min_sample_size2: int = 3,
):
    """
    Compute the spearman correlation for two independent samples.

    Parameters
    ----------
    data : pd.DataFrame
        Data in "record" format. Each row represents a record. Indices represent groups.
    groupby : str
        The name of the column that contains the group information
    group1 : str
        The name of the first group
    group2 : str
        The name of the second group
    min_sample_size1 : int
        The minimal sample size of the first group
    min_sample_size2 : int
        The minimal sample size of the second group
    value : str
        The name of the column that contains the value information

    Returns
    -------
    statistic
        The spearman correlation of the test
    pval
        The pval of the test
    mean
        The mean of each group
    std
        The standard deviation of each group
    median
        The median of each group
    sample_size
        The sample size of each group
    """

    # get the data
    data1 = data.query(f"`{groupby}` == '{group1}'")[value].dropna()
    data2 = data.query(f"`{groupby}` == '{group2}'")[value].dropna()

    # get the sample size
    sample_size = [len(data1), len(data2)]

    # calculate the mean
    mean = [np.mean(data1), np.mean(data2)]
    std = [np.std(data1), np.std(data2)]
    median = [np.median(data1), np.median(data2)]

    # perform the test
    if len(data1) >= min_sample_size1 and len(data2) >= min_sample_size2:
        statistic, pval = spearmanr(data1, data2)
    else:
        statistic, pval = np.nan, np.nan

    # return the result
    return statistic, pval, mean, std, median, sample_size


def H_test(data):
    """
    Compute the Kruskal-Wallis H-test for independent samples. It is a non-parametric version of ANOVA.
    See details, https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kruskal.html

    Parameters
    ----------
    data : pd.DataFrame
        Data in "record" format. Each row represents a record. Indices represent groups.

    Returns
    -------
    statistic
        The H statistic
    pval
        The pval for the test
    """

    tmp = data.sort_index()
    tmp = tmp[~tmp.index.isna()]  # remove group of nan

    if len(tmp.index.unique()) > 1:
        mat = []
        for idx in tmp.index.unique():
            mat.append(tmp.loc[idx].values.flatten())
        statistic, pval = kruskal(*mat)
    else:
        statistic, pval = np.nan, np.nan

    return statistic, pval


def f_oneway_test(data):
    """
    Perform one-way ANOVA.
    See details, https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.f_oneway.html

    Parameters
    ----------
    data : pd.DataFrame
        Data in "record" format. Each row represents a record. Indices represent groups.

    Returns
    -------
    statistic
        The F statistic of the test
    pval
        The pval for the test
    """
    tmp = data.sort_index()
    tmp = tmp[~tmp.index.isna()]  # remove group of nan

    if len(tmp.index.unique()) > 1:
        mat = []
        for idx in tmp.index.unique():
            mat.append(tmp.loc[idx].values.flatten())
        statistic, pval = f_oneway(*mat)
    else:
        statistic, pval = np.nan, np.nan

    return statistic, pval
