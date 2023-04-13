import matplotlib as mpl
import matplotlib.pyplot as plt

class Plotter:

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
