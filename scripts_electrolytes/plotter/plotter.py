import matplotlib as mpl
import matplotlib.pyplot as plt
import os

plt.rc('text', usetex=True)

class Plotter:

    '''
        Kwargs options:
            kwargs with default values:
                figsize: figure size, in inches
                linewidth: linewidth for Line2D objects
                labelsize: axes labels fontsize
                figname: filename for the figure. 
                savefig: Boolean, should the figure be saved or not.
                showfig: Boolean, display figure or not

            kwargs without default values (ignored if not specified):
                title: figure title
    '''

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

    def set_line2d_params(self, defname, **kwargs):

        self.set_linewidth(kwargs.get('linewidth', 2.0))
        self.set_labelsize(kwargs.get('labelsize', 16))
        self.set_figname(kwargs.get('figname', defname))
        self.set_savefig(kwargs.get('savefig', True))
        self.set_showfig(kwargs.get('showfig', True))
        self.set_marker(kwargs.get('marker', 'o'))
        self.set_linecolor(kwargs.get('linecolor', False))
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

    def set_savefig(self, savefig):
        self.savefig = savefig

    def set_title(self, title):
        self.title = title

    def set_showfig(self, showfig):
        self.showfig = showfig

    def set_marker(self, marker):
        self.marker = marker

    def set_linecolor(self, col):
        self.linecolor = col

    def add_title(self):
        self.ax.set_title(self.title, fontsize=self.labelsize)

    def save_figure(self):

        if self.savefig:
            try:
                os.mkdir('figures')
            except OSError:
                pass
            pathname = os.path.join('figures', self.figname)
            plt.savefig(pathname, bbox_inches='tight')

    def show_figure(self):
        if self.showfig:
            plt.show()
