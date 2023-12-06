import os
from ..plotter.neb_plotter import NebPathPlotter
from ..plotter.colorpalettes import bright
from ase.io import write as ase_write
import warnings

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


class NebTraj:

    def __init__(self, nreplica):

        self.nreplica = nreplica


    def write_neb_trajectory(self, out_fname='nebtraj.xyz'):

        self.check_output_format(out_fname)
        ase_write(filename=self.out_fname, images=self.trajectory, append=True)


    def check_output_format(self, fname):
        fmt = os.path.splitext(fname)[1]
        # In multi-configuration write mode, ASE accepts cif and xyz
        # Although I removed cif since Ovito does not accept it in multi-config mode
        # I could also add .traj, but Ovito requires the paid version for this type of file...
        
        # For now, force output in xyz format
        if fmt != '.xyz':
            msg = '''The only output format supported by ASE and by Ovito in multiconfiguration mode is .xyz, so I will output the trajectory in this format instead of {}.'''.format(fmt)
            warnings.warn(msg)
            self.out_fname = os.path.splitext(fname)[0]+'.xyz'
        else:
            self.out_fname = fname

        # fmt_is_ok = ['.xyz']

        # if fmt not in fmt_is_ok:
        #     raise ValueError('output format should be in {} but I got {}'.format(fmt_is_ok, fmt))
