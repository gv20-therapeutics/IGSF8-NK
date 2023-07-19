from typing import Any, Dict, Optional
import matplotlib.pylab as pylab
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from OncoMarkAI.utils import get_logger


# Figure setting function. A default figure setting is provided. You can also update the default setting by providing a dictionary of parameters.
def fig_setting(
    params: Optional[Dict[str, Any]] = None,
):
    """
    Set figure parameters
    """
    # default params
    default_params = {
        "figure.facecolor": "w",
        "figure.figsize": (4, 3),
        "figure.autolayout": True,
        "figure.dpi": 300,
        "axes.titlesize": 6,
        "axes.linewidth": 0.5,
        "axes.grid": False,
        "xtick.labelsize": 5,
        "ytick.labelsize": 5,
        "legend.fontsize": 4,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "font.size": 6,
        "pdf.fonttype": 42,
        "mathtext.default": "regular",
    }

    # update default params
    if params is not None:
        default_params.update(params)
    pylab.rcParams.update(default_params)

    # seaborn setting
    # sns.set_palette("muted")
    # sns.set_style("ticks")
    # sns.despine(offset=10, trim=True)

    # logger
    plot_logger = get_logger("plotting")
    plot_logger.info("Figure setting updated")


def pval2star(pval, keep_pval=True):
    """
    Generate a string of pval stars

    Parameters
    ----------
    pval : float
        pval
    keep_pval : bool
        whether to keep the pval in the string
    Returns
    -------
    star+pval : str
    """

    if pval <= 0.0001:
        star = "****"
    elif pval <= 0.001:
        star = " ***"
    elif pval <= 0.01:
        star = " **"
    elif pval <= 0.05:
        star = "  *"
    else:
        star = "ns"

    if keep_pval:
        star += f" p={pval:.1e}"
    else:
        star = star.replace("ns", " ")

    return star


class cmap:
    """
    Create a colorblind palette of publication quality
    """

    def __init__(self) -> None:
        colorblind_palette = {"Orange": "#E69F00", "Skyblue": "#56B4E9", "Green": "#009E73", "Reddish_purple": "#CC79A7", "Yellow": "#F0E442", "Blue": "#0072B2", "Vermilion": "#D55E00", "Gray": "#999999", "Black": "#000000", "White": "#FFFFFF"}
        self.colorblind_palette = ListedColormap(list(colorblind_palette.values()), name="colorblind_palette")

        npg_palette = {"Orangered": "#E64B35FF", "Deepskyblue": "#4DBBD5FF", "Seagreen": "#00A087FF", "Steelblue": "#3C5488FF", 
                       "Coral": "#F39B7FFF", "Lightsteelblue": "#8491B4FF", "Turquoise": "#91D1C2FF", "Red": "#DC0000FF", "Saddlebrown": "#7E6148FF",
                       "darkgoldenrod": "#B09C85FF", "White": "#FFFFFF", "Black": "#000000"}  
        self.npg_palette = ListedColormap(list(npg_palette.values()), name="npg_palette")

    def __getitem__(self, key):
        return self.__dict__[key]
    