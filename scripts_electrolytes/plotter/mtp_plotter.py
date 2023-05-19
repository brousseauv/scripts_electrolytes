import matplotlib.pyplot as plt
from ..interfaces.mtp_interface import read_errors
from .plotter import Plotter

''' Base class for plotting validation and test errors on MTP models.

    Input:
        ax: matplotlib.axes.AxesSubplot instance 
        Default is None, in which case a new figure/ax is created

    Output: matplotib.axes.AxesSubplot instance containing the data extracted from given MTP output file

'''
# define functions that return a matplotlib ax instance, then is can freely be customized by the user
# Just as I did in ThermalExpansion.teplotter

class MtpPlotter(Plotter):
    
    def __init__(self, ax=None, **kwargs):

        self._init_figure(ax=ax, **kwargs)


def check_field(field):

    field_list = ['energy', 'forces', 'stresses', 'all']
    if field not in field_list:
        raise Exception('Field value should one of {} but I got {}'.format(field_list, field))

# I can make different ploting functions depending on what I'm looking into
# i.e. errors, convergence/active learning stuff... 
# I leave this open for future needs
# have a variable to pass errors ans errorbars (for multiple models, if needed

def plot_errors(
        ax=None,
        field='all',
        data=None,
        runtype='test',
        label_size=16,

        **kwargs):

    # FIX ME: right now this is sort of useless unless I can track error evolution through an active learning scheme
    # since there are no epochs here, there is only one error report per passive run. 
    check_field(field)

    # set n_subplots in case axes are not predefined, need to be able to handle this when creating figure from scratch

    # if ax return ax else figure.show (add function in Plotter)
