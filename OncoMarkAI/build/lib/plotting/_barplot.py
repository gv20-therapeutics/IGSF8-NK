from typing import Optional, Union, List
from itertools import combinations
import matplotlib
import seaborn as sns
import anndata as ad

from ..tools import _stats as stats
from ._setting import pval2star


def barplot(
    adata: ad.AnnData,
    x: str,
    y: str,
    order: Union[List[str], str] = None,
    groupby: Optional[str] = None,
    palette: Optional[List[str]] = None,
    color: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    **kwargs,
) -> matplotlib.axes.Axes:
    """Boxplot to show the distribution of a variable.

    Parameters
    ----------
    adata : ad.AnnData
        TCGA data saved in AnnData format.
    x : str
        The variable to be plotted on the x axis, e.g. "tumor types".
    y : str
        The variable to be plotted on the y axis, e.g. gene expression and "MSI sensor score".
    order : Optional[List[str], str], optional
        The order of x variable. Can be a list. If 'decending', the order will be descending. By default None. If None, the order will be the unique values in the x axis.
    groupby : Optional[str], optional
        The variable to groupby, by default None
    jitter : bool, optional
        Whether to plot strippplot, by default False
    palette : Optional[List[str]], optional
        Color palette, by default None
    color : Optional[str], optional
        Color of the boxplot, by default None
    xlabel : Optional[str], optional
        Xlable to show on the xaxis, by default None
    ylabel : Optional[str], optional
        Ylabel to show on the yaxis, by default None
    title : Optional[str], optional
        The plot title, by default None
    width : float, optional
        The width of the box, by default 0.65
    add_stats : bool, optional
        Whether to add the pvalue, by default False

    Returns
    -------
    matplotlib.axes.Axes
        The plot object.

    """
    # check if the x axis is in the adata.obs
    if x not in adata.obs.columns:
        raise ValueError(f"{x} is not in adata.obs.columns")

    # check if the y axis is in the adata.obs or adata.var
    if y not in adata.obs.columns and y not in adata.var_names:
        raise ValueError(f"{y} is not in adata.obs.columns or adata.var.columns")

    # check if the y is gene expression
    if y not in adata.obs.columns and y in adata.var_names:
        adata.obs[y] = adata[:, y].to_df()

    # check if the groupby is provided
    if groupby is None:
        data = adata.obs[[x, y]]
    else:
        if groupby not in adata.obs.columns:
            raise ValueError(f"{groupby} is not in adata.obs.columns")
        data = adata.obs[[x, y, groupby]]

    # check if the order is provided
    if order is None:
        order = list(adata.obs[x].unique())
    elif not isinstance(order, list):
        if isinstance(order, str) and order == "decending":
            order = (
                data.groupby(x).mean().sort_values(by=y, ascending=False).index.tolist()
            )

    # check if the palette is provided
    # if palette is None:
    #     palette = sns.color_palette(n_colors=len(order))

    sns.set_style("ticks")
    ax = sns.catplot(
        data=data,
        x=x,
        order=order,
        y=y,
        hue=groupby,
        kind="bar",
        palette=palette,
        color=color,
        legend=False,
        **kwargs,
    )

    sns.despine(top=False, right=False, left=False, bottom=False)

    if xlabel is not None:
        ax.set(xlabel=xlabel)
    if ylabel is not None:
        ax.set(ylabel=ylabel)
    if title is not None:
        ax.set(title=title)

    return ax
