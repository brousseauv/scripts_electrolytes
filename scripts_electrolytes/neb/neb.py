import os
from ..plotter.neb_plotter import NebPathPlotter
from ..plotter.colorpalettes import bright

''' Basic class that calls the plotter '''
# this is a different class in case I add Abinit output possibilities

class NebData:

    def __init__(self, fname, rootname):

        self.fname = fname
        self.rootname = rootname


    def plot_neb_energy(self, plot_verbose=True, **kwargs):

        myplot = NebPathPlotter(**kwargs)
        myplot.set_line2d_params(defname='neb.png', **kwargs)

        if myplot.linecolor:
            col = myplot.linecolor
        else:
            col = bright['blue']

        myplot.ax.plot(self.reaction_coordinate, self.potential_energy, color=col, marker=myplot.marker,
                       linestyle=myplot.linestyle, linewidth=myplot.linewidth)

        if plot_verbose:
            myplot.ax.text(0.05, 0.90, r'E$_A$={:.3f} eV'.format(self.forward_barrier), fontsize=myplot.labelsize+1, transform=myplot.ax.transAxes)
        try:
            myplot.add_title()
        except:
            pass

        myplot.set_labels()
        myplot.save_figure()
        myplot.show_figure()


    def print_barriers(self):

        print('Forward (backward) energy barrier: {:.4f} ({:.4f}) eV'.format(self.forward_barrier, self.backward_barrier))
