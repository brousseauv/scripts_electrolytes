import matplotlib as mpl
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)

class Plotter:

    def __init__(self, ax=None, **kwargs):

        self._init_figure(ax=ax, **kwargs)

    def _init_figure(self, ax=None, **kwargs):

        if ax:
            self.ax = ax
            self.fig = ax.get_figure()
        else:
            self.fig, self.ax = plt.subplots(1, 1)
            self.set_size(*kwargs.get('figsize', (6, 6)))
            # FIX ME: test the code without providing ax in the input

    def set_size(self, w, h):
        ''' Set figure size'''
        self.fig.set_size_inches(w, h)

    def set_line2d_params(self,**kwargs):

        self.set_linewidth(kwargs.get('linewidth', 2.0))
        self.set_labelsize(kwargs.get('labelsize', 16))
        self.set_figname(kwargs.get('savefig', 'diffusion.png'))
        self.set_showfig(kwargs.get('showfig', True))
        try:
            self.set_title(kwargs.get('title'))
        except:
            pass

    def set_linewidth(self, lw):
        self.linewidth = lw

    def set_labelsize(self, ls):
        self.labelsize = ls

    def set_figname(self, name):
        self.figname = name

    def set_title(self, title):
        self.title = title

    def set_labels(self):
        self.ax.set_xlabel(r'time (ps)', fontsize=self.labelsize)
        self.ax.set_ylabel(r'MSD (\AA$^2$)', fontsize=self.labelsize)

    def set_showfig(self, showfig):
        self.showfig = showfig

    def set_title(self):
        self.ax.set_title(self.title)

    def save_figure(self):
        plt.savefig(self.figname)

    def show_figure(self):
        if self.showfig:
            plt.show()
