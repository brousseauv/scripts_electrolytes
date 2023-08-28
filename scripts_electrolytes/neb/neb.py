import os
from ..plotter.neb_plotter import NebPathPlotter
from ..plotter.colorpalettes import bright

''' Basic class that calls the plotter '''
# this is a different class in case I add Abinit output possibilities

class NebData:

    def __init__(self, fname, rootname):
        try:
            os.path.exists(fname)
        except:
            raise NameError('File {} does not exists. Please provide correct path to the "log.lammps" output file containing NEB results.'.format(fname))

        self.fname = fname
        self.rootname = rootname


    def plot_neb_energy(self, verbose=True, **kwargs):

        myplot = NebPathPlotter(**kwargs)
        myplot.set_line2d_params(defname='neb.png', **kwargs)

        if myplot.linecolor:
            col = myplot.linecolor
        else:
            col = bright['blue']

        myplot.ax.plot(self.reaction_coordinate, self.potential_energy, color=col, marker=myplot.marker,
                       linestyle='solid', linewidth=myplot.linewidth)

        if verbose:
            myplot.ax.text(0.05, 0.90, r'E$_A$={:.3f} eV'.format(self.forward_barrier), fontsize=myplot.labelsize+1, transform=myplot.ax.transAxes)
        try:
            myplot.add_title()
        except:
            pass

        myplot.set_labels()
        myplot.save_figure()
        myplot.show_figure()
