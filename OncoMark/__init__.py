import sys
from importlib.metadata import version

package_name = "OncoMark"
__version__ = version(package_name)

from . import tools as tl
from . import plotting as pl
from . import data

sys.modules.update({f"{__name__}.{m}": globals()[m] for m in ["tl", "pl"]})
# clean up namespace
for _ in ["tools", "plotting", "sys", "utils", "version"]:
    del globals()[_]
