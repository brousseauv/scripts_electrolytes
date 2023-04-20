from ..interfaces.lammps_interface import extract_thermo, create_thermo_dataframe, read_msd
from .msd import MsdData

class LammpsMsdData(MsdData):

    def __init__(self, fname, filetype):

        if filetype not in ['thermo', 'dump', 'dump-netcdf']:
            raise Exception('''filetype should be one of the following: "thermo", "dump",
                               "dump-netcdf", but I got {}'''.format(filetype))

        self.filetype = filetype
        super(LammpsMsdData, self).__init__(fname)

    
    def compute_diffusion_coefficient(self, thermo_fname=None, plot=False, **kwargs):

        if self.filetype == 'thermo':
            if not thermo_fname:
                thermo_fname = extract_thermo(self.fname)
            data = create_thermo_dataframe(thermo_fname)
            self.time, self.msd = read_msd(data)

        elif self.filetype == 'dump':
            raise NotImplementedError('"dump" filetype not yet implemented')

        elif self.filetype == 'dump-netcdf':
            raise NotImplementedError('"dump-netcdf" filetype not yet implemented')

        diffusion = self.extract_diffusion_coefficient()
        print('Diffusion coefficient: {:.3e} cm^2/s'.format(diffusion))

        if plot:
            self.plot_diffusion_coefficient(**kwargs)
