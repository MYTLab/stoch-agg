import os
import errno
import subprocess as sp


class KineticMC:
    """  """
    para_template = """nsample {nsample}
nstep 10000000
nmonomer {nmonomer}
nc {nc}
kn {kn}
kp {kp}
kf {kf}
km {km}
k2 {k2}
n2 {n2}
Nt {Nt}
Ndatapoint {Ndatapoint}
p0 {p0}
{data_path}"""

    def __init__(self, data_folder=None, exe_path=None, exe_name='AG', **para):
        self.exe_path = os.path.join(os.getcwd(), '..', 'kMC') if exe_path is None else exe_path
        self.exe_name = exe_name
        self.data_path = os.path.join('..', 'kMC_data', data_folder, '') if data_folder is None else data_folder
        self.para = para

    def run(self, cluster=None, verbose=False, overwrite=False):
        """ Main method. Wrapper for creating parameters.ini and running simulation.
            Set check to False when running in non-interactive mode, such as on a cluster.
            cluster = 'qsub test.sh'
         """
        # check exe exist
        if not self.check_path_exit(os.path.join(self.exe_path, self.exe_name)):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), self.exe_name)
        # check data folder exist
        if self.check_path_exit(self.data_path):
            msg = 'Data folder exists. Stop to prevent overwriting of data existed.'
            assert overwrite, msg
            proc = sp.Popen(['rm -r ' + self.data_path], shell=True, stdout=sp.PIPE)
            proc.communicate()
        os.mkdir(self.data_path)
        print(f'Folder {self.data_path} created.')

        self.write_paraini()
        if cluster is None:
            pout, perr = self.execute()
            st = 'finished'
        else:
            pout, perr = self.submit(*cluster.split(' '))
            st = 'submitted'
        if perr:
            print('Simulation failed.')
            print(perr)
        if verbose:
            print(pout)
        
        print('Simulation '+ st)

    # TODO:
    def write_submit_script(self):
        pass

    def write_paraini(self):
        """ Create parameters.ini in the same dir as the executable """
        p = {'data_path': self.data_path}
        p.update(self.para)
        file_name = os.path.join(self.exe_path, 'parameters.ini')
        with open(file_name, 'w+') as f:
            f.write(self.para_template.format(**p))
            print(f'{file_name} created.')

    def execute(self):
        """ Run the kmc executable. Return std out and err. """
        dir_n_exe = self.exe_path + ' && ./AG parameters.ini'
        proc = sp.Popen(['cd ' + dir_n_exe], shell=True, stdout=sp.PIPE)
        pout, perr = proc.communicate()
        return pout, perr
    
    def submit(self, submit_cmd, submit_script):
        """ Submit the task to a cluster.  """
        dir_n_sub = self.exe_path + f' && {submit_cmd} {submit_script}'
        proc = sp.Popen(['cd ' + dir_n_sub], shell=True, stdout=sp.PIPE)
        pout, perr = proc.communicate()
        return pout, perr

    @staticmethod
    def check_path_exit(path):
        return os.path.exists(path)



if __name__ == '__main__':
    test_para = dict(nsample=5, p0=0, nc=2, n2=0, kn=1e-7, kp=1e1,
                     km=1e-6, k2=0, kf=1e-8, nmonomer=100, Nt=5, Ndatapoint=100)
#    test = KineticMC('test', **test_para)
#    test.run()
