# Stochastic Protein Aggregation
import numpy as np
import Stochastic_model
from scipy.integrate import odeint
import pickle
import pandas as pd
import os
import matplotlib.pyplot as plt
from kineticMonteCarlo import KineticMC

# Author: Ace Shen, 2019-10-13

class StochasticAgg:
    """"""
    valid_parameters = ['nc', 'n2', 'kn', 'kp', 'km', 'kf', 'k2', 'nmonomer', 'p0', 'Nt', 'Ndatapoint']

    def __init__(self, verbose=False):
        self.para = {i: None for i in self.valid_parameters}
        self.z, self.z0, self.dzdt  = None, None, None
        self.kmc = None
        self.stochmodel = False
        self.func_name = None
        self.kmc_touched = False
        self.verbose = verbose
        self.closeorder = None

        self.kmc_pmean, self.kmc_pstd, self.kmc_mmean, self.kmc_mstd = None, None, None, None
        self.kmc_p, self.kmc_m = None, None
        self.ana_pmean, self.ana_mmean, self.ana_pstd, self.ana_mstd = None, None, None, None
        self.t = None
    
        self.kmc_dpdt, self.kmc_dmdt = None, None
        self.ana_dpdt, self.ana_dmdt = None, None

    def __repr__(self):
        out0 = ('kMC: ' + str(self.kmc)) if self.kmc_touched else 'kMC: None'
        out1 = 'closeorder: ' + str(self.closeorder) if self.closeorder is not None else 'closeorder: None'
        out2 = 'para: ' + str(self.para)
        kmc = '\navailable kmc data: kmc_pmean, kmc_pstd, kmc_mmean, kmc_mstd, kmc_dpdt, kmc_dmdt' if self.kmc_touched else ''
        ana = '\navailable ana data: ana_pmean, ana_pstd, ana_mmean, ana_mstd, ana_dpdt, ana_dmdt' if self.stochmodel else ''
        return out0 + '\n' + out1 + '\n' + out2 + kmc + ana
    
    def _printout(self, s, end='\n', v=False):
        if v or self.verbose:
            print(s, end=end)
    
    # ====================================================================================
    # ============== Parameters ========================================================
    # ====================================================================================

    def set_para(self, from_kmc=None, closeorder=None, **para):
        """ Set parameters
            If from_kmc is True, get parameters from the simulation output 'para.txt'.
            Otherwise pass a dict with parameters.
            """
        if from_kmc:
            self.kmc = from_kmc
            self._get_para_from_kmc(from_kmc)
            self.kmc_touched = False
        else:
            input = [i for i in list(para.keys()) if i not in self.valid_parameters]
            assert not input, f'{input} are not valid parameters.'
            self.para.update(para)
        # time axis
        if (self.para['Nt'] is not None) & (self.para['Ndatapoint'] is not None):
            self._set_time(self.para['Nt'], self.para['Ndatapoint'])
        # initial m
        self.para['m0'] = np.ceil(self.para['nmonomer'] - self.para['nc'] * self.para['p0'])
        if self.para['p0'] == 0:
            self.para['p0'] += 1e-6
        self._set_stochmodel(closeorder)
        self.stochmodel = False  # always set this to False after updating parameters

    def _set_stochmodel(self, closeorder):
        self.closeorder = closeorder
        if closeorder is not None:
            if closeorder == 0:
                self.func_name = 'deterministic'
            else:
                self.func_name = 'stochastic_close' + str(self.closeorder)\
                                + '_nc' + str(self.para['nc'])\
                                + '_n2' + str(self.para['n2'])
            self.model = getattr(Stochastic_model, self.func_name)
            self._ordered_para = (self.para['kn'], self.para['kp'],
                                  self.para['km'], self.para['kf'],
                                  self.para['nmonomer'], self.para['k2'],
                                  self.para['nc'], self.para['n2'])
        else:
            print('Closeorder is not set. Using ONLY for parsing kMC data.')
            print('To set closeorder later, use set_para(closeorder=c)')

    def _get_para_from_kmc(self, kmc_path):
        keys_int = ['nsample', 'nmonomer', 'nc', 'n2', 'Ndatapoint', 'p0']
        keys_f = ['kf', 'km', 'k2', 'Nt']
        special_keys = ['kn [1/s]', 'kp [1/s]']
        try:
            with open(os.path.join(kmc_path, 'para.txt'), 'r') as f:
                for line in f:
                    ll = line.strip().split(" = ")
                    if ll[0] in keys_int:
                        self.para[ll[0]] = int(ll[1])
                    elif ll[0] in keys_f:
                        self.para[ll[0]] = float(ll[1])
                    elif ll[0] in special_keys:
                        self.para[ll[0][:2]] = float(ll[1])
                    else:
                        pass
        except FileNotFoundError as err:
            print(err)

    def _set_time(self, end, n, start=0, recursive=None):
        # set time axis for solution
        # if used in recursive mode, return a time slice for the use of odeint
        self.t = np.linspace(start, end, n+1, endpoint=True)
        if recursive is not None:
            return self.t[-recursive-1:]
    
    # ====================================================================================
    # ============== kMC data ========================================================
    # ====================================================================================

    def run_cppkmc(self, data_folder):
        """ Run the kMC simulation through the cpp program """
        simulation = KineticMC(data_folder, **self.para)
        simulation.run()
        # set kmc data path here...

    def runkmc(self):
        """ Parse kMC data """
        assert self.kmc, 'Need kmcdata path in StochasticAgg.set_para()'
        self._printout(f'Processing kmc data from {self.kmc} ', end='...')
        self.kmc_touched = True
        #self._load_kmc_maxl()
        self._load_kmc_time()
        self._load_kmcdata()
        self._process_kmcdata()
        self._printout(' Done.')
    
    def _load_kmc_maxl(self):
        # Read in maxl, get the number of how many files to read for each sample
        try:
            with open(os.path.join(self.kmc, 'maxl.txt'), 'r') as f:
                self._maxl = np.loadtxt(f, dtype=int).reshape(-1)  # make sure its an array of sequence
        except FileNotFoundError as err:
            print(err)

    def _load_kmc_time(self):
        # Read in time
        try:
            with open(os.path.join(self.kmc, 'time.dat'), 'r') as f:
                self.kmc_t = np.loadtxt(f).reshape(-1)  # make sure its an array of sequence
        except FileNotFoundError as err:
            print(err)

    def _load_kmcdata(self):
        nsample = self.para['nsample']
        ndatapoint = self.para['Ndatapoint']
        nmonomer = self.para['nmonomer']
        # get p(t) and m(t) from each sample
        self.kmc_p = np.zeros((nsample, ndatapoint+1))
        self.kmc_m = np.zeros((nsample, ndatapoint+1))
        # data filename
        filename, filetype = 'adata_', '.dat'
        # loop over each sample, get all the l in multiple files for each sample
        for isample in range(0, nsample):
            # read data
            fn = os.path.join(self.kmc, filename + str(isample+1) + filetype)
            try:
                with open(fn, 'r') as f:
                    data = np.loadtxt(f)
            except FileNotFoundError as err:
                print(err)
            # First moment
            self.kmc_p[isample] = np.sum(data[1:, :], axis=0)
            # Second moment = total number of monomers in aggregates
            self.kmc_m[isample] = nmonomer - data[0]

    def _process_kmcdata(self):
        msg1 = f'Length of p in time is {self.kmc_p.shape[1]} but t is {len(self.t)}'
        msg2 = f'Length of m in time is {self.kmc_m.shape[1]} but t is {len(self.t)}'
        assert self.kmc_p.shape[1] == len(self.t), msg1
        assert self.kmc_m.shape[1] == len(self.t), msg2
        self.kmc_pmean = self.kmc_p.mean(axis=0)
        self.kmc_pstd = self.kmc_p.std(axis=0)
        self.kmc_mmean = self.kmc_m.mean(axis=0)
        self.kmc_mstd = self.kmc_m.std(axis=0)
    
        self.kmc_dpdt = np.gradient(self.kmc_pmean, self.t)
        self.kmc_dmdt = np.gradient(self.kmc_mmean, self.t)
    
    # ====================================================================================
    # ============== Stochastic Model ========================================================
    # ====================================================================================
    
    def run(self, recursive=None):
        """ Solve the stochastic model
            Need to set closeorder before running. If kmc is set in self.set_para, parse kMC data.
            """
        self._set_initial_condition(recursive=recursive)
        self._solve_eq(recursive=recursive)
        if recursive is None:
            if (self.kmc is not None) & (self.kmc_touched is False):
                self.runkmc()
            self._parse_solution()
            self._get_slope()
            self.stochmodel = True
    
    def _set_initial_condition(self, recursive=False):
        err_msg = 'Need to set closeorder for a stochastic model => StochasticAgg(closeorder=c)'
        assert self.closeorder is not None, err_msg
        if recursive is not None:
            self.z0 = self.z[-1, :]
        else:
            if self.closeorder == 0:
                self.z0 = [self.para['p0'], self.para['m0']]
            else:
                self.z0 = [pow(self.para['p0'], iplusj - j)*pow(self.para['m0'], j) for iplusj in range(1, self.closeorder+1) for j in range(iplusj+1)]

    def _solve_eq(self, recursive=None):
        # solving the stochastic model
        if recursive is None:
            self.z = odeint(self.model, self.z0, self.t, args=(self._ordered_para,))
        else:
            n_recur, nd, t_half = recursive[0], recursive[1], recursive[2]
            ode_t = self._set_time(n_recur*t_half, n_recur*nd, recursive=nd)
            self.z = np.vstack(( self.z[:-1, :], odeint(self.model, self.z0, ode_t, args=(self._ordered_para,)) ))

    def _parse_solution(self):
        # z: [p1m0, p0m1, p2m0, p1m1, p0m2]
        self.ana_pmean = self.z[:, 0]
        self.ana_mmean = self.para['nmonomer'] - self.z[:, 1]
        if self.closeorder != 0:
            ana_pvar = self.z[:, 2] - self.ana_pmean*self.ana_pmean
            self.ana_pstd = pow(ana_pvar, 0.5)
            ana_mvar = self.z[:, 4] - self.z[:, 1]*self.z[:, 1]
            self.ana_mstd = pow(ana_mvar, 0.5)

    def _get_slope(self):
        # all models are not explicitly time-dependent
        self.dzdt = np.empty_like(self.z)
        for i in range(len(self.dzdt)):
            self.dzdt[i] = self.model(self.z[i], 0, self._ordered_para)
        
        self.ana_dpdt = self.dzdt[:, 0]
        self.ana_dmdt = -self.dzdt[:, 1]

    # ====================================================================================
    # ============== Steady solution ========================================================
    # ====================================================================================

    def solve_steady(self, m_ths=0.1, nd=100, iter_nochange=10, max_iter=500):
        """ Find steady state of m. Stop either m <= m_ths or m doesn't change for iter_nochange iterations.
            Nd: The number of time steps in each iteration.
            """
        # Using t_half formula from:
        # Mean-field master equation formalism for biofilament growth
        # Thomas C. T. Michaels and Tuomas P. J. Knowles
        # American Journal of Physics 82, 476 (2014)
        estimate_thalf = (1/(2*self.para['kn'] * self.para['kp'] * self.para['m0']**self.para['nc'])**0.5)
        self._printout(f'The estimated t_half is {estimate_thalf} s.')
        self._set_time(estimate_thalf, nd)
        check_nochange = np.zeros(iter_nochange)
        # first solve
        self.run()
        count = 1
        m = self.z[-1, 1]/self.para['nmonomer']  # free monomer
        check_nochange[-1] = m
        while m > m_ths:
            count += 1
            self.run(recursive=(count, nd, estimate_thalf))
            m = self.z[-1, 1]/self.para['nmonomer']
            check_nochange = pd.Series(check_nochange).shift(-1).to_numpy()
            check_nochange[-1] = m
            round = np.around(check_nochange, decimals=4)
            if all(round - round[0] == 0):
                self._printout(f'Steady state found after {count} iterations.')
                break
            if count > max_iter:
                self._printout(str(self.para), v=True)
                self._printout(f'Maximum iteration exceeded. max_iter = {max_iter}.', v=True)
                break
                
        self.para['Nt'], self.para['Ndatapoint'] = self.t[-1], count*nd
        self._printout(f'Finished. Final free monomer percentage: {m:.3f}', end='; ')
        self._printout(f'Stop time: {self.t[-1]} s')
        self._parse_solution()
        self._get_slope()
        self.stochmodel = True

    # ====================================================================================
    # ============== Pickle data ========================================================
    # ====================================================================================

    def to_pickle(self, filename, path=None, deep=False):
        """ Save the analysis into a pkl file. """
        # There is also no guarantee of compatibility between different versions of Python
        # because not every Python data structure can be serialized by the module.

        if not path:
            path = os.getcwd()
        r = os.path.join(path, filename)
        p = {}
        if self.kmc:
            self._pickle_kmc(p, deep=deep)
        if self.stochmodel:
            self._pickle_model(p)
        with open(r, 'wb') as dumpfile:
            pickle.dump(p, dumpfile)
        self._printout('Data pickled to ' + r)

    def from_pickle(self, filename, path=None):
        """ Read previous analysis from file. """
        if not path:
            path = os.getcwd()
        with open(os.path.join(path, filename), 'rb') as readfile:
            p = pickle.load(readfile)
        self.para = p['parameters']
        if 'kmc_datapath' in list(p.keys()):
            self.kmc_touched = True
            self.kmc = p['kmc_datapath']
            self.para = p['parameters']
            self.kmc_pmean = p['kmc_pmean']
            self.kmc_pstd = p['kmc_pstd']
            self.kmc_mmean = p['kmc_mmean']
            self.kmc_mstd = p['kmc_mstd']
            self.kmc_dpdt = p['kmc_dpdt']
            self.kmc_dmdt = p['kmc_dmdt']
            self.kmc_t = p['kmc_t']
            if 'kmc_p' in list(p.keys()):
                self.kmc_p = p['kmc_p']
                self.kmc_m = p['kmc_m']
        
        if 'func_name' in list(p.keys()):
            self.stochmodel = True
            self.func_name = p['func_name']
            self.ana_pmean = p['ana_pmean']
            self.ana_mmean = p['ana_mmean']
            self.ana_dpdt = p['ana_dpdt']
            self.ana_dmdt = p['ana_dmdt']
            self.closeorder = p['closeorder']
            self.t = p['t']
            if self.closeorder != 0:
                self.ana_pstd = p['ana_pstd']
                self.ana_mstd = p['ana_mstd']
            self._set_stochmodel(self.closeorder)

    def _pickle_kmc(self, p, deep):
        p['kmc_datapath'] = self.kmc
        p['parameters'] = self.para
        p['kmc_pmean'] = self.kmc_pmean
        p['kmc_pstd'] = self.kmc_pstd
        p['kmc_mmean'] = self.kmc_mmean
        p['kmc_mstd'] = self.kmc_mstd
        p['kmc_dpdt'] = self.kmc_dpdt
        p['kmc_dmdt'] = self.kmc_dmdt
        p['kmc_t'] = self.kmc_t
        if deep:
            p['kmc_p'] = self.kmc_p
            p['kmc_m'] = self.kmc_m


    def _pickle_model(self, p):
        p['func_name'] = self.func_name
        p['parameters'] = self.para
        p['ana_pmean'] = self.ana_pmean
        p['ana_mmean'] = self.ana_mmean
        p['ana_dpdt'] = self.ana_dpdt
        p['ana_dmdt'] = self.ana_dmdt
        p['closeorder'] = self.closeorder
        p['t'] = self.t
        if self.closeorder != 0:
            p['ana_pstd'] = self.ana_pstd
            p['ana_mstd'] = self.ana_mstd


    # ====================================================================================
    # ============== Plot ========================================================
    # ====================================================================================

    def plot(self):
        """ Plot all results if they exist.
            Return two figure objects for P(t) and M(t).
            """
        figp, axp = plt.subplots(figsize=(8, 6))
        figm, axm = plt.subplots(figsize=(8, 6))
        title = ''
        if self.stochmodel:
            self.plot_ana([axp, axm])
            title += 'closeorder ' + str(self.closeorder)
        if self.kmc:
            self.plot_kmc([axp, axm])
            title += '  ' + str(self.para['nsample']) + ' samples'
        if (not self.kmc) and (not self.stochmodel):
            print('No analysis performed yet.')
            return
        axp.set_title( title, fontsize=20)
        axm.set_title( title, fontsize=20)
        return figp, figm
    
    def plot_ana(self, ax=None):
        """ Plot analytical result only.
            Return two figure objects if parameter \'ax\' is not provided.
            """
        if ax:
            axp, axm = ax[0], ax[1]
            figp, figm = None, None
        else:
            figp, axp = plt.subplots(figsize=(8, 6))
            figm, axm = plt.subplots(figsize=(8, 6))
        t = self.t
        # plt.fill_between(tode, ana_mmean+ana_mstd, ana_mmean-ana_mstd, alpha=0.5, label = r'$1\sigma$-Ana')

        axp.plot(t, self.ana_pmean, '-', color='blue', lw=3, label='model')
        if self.closeorder != 0:
            axp.fill_between(t, self.ana_pmean + self.ana_pstd,
                             self.ana_pmean - self.ana_pstd, color='dodgerblue', alpha=0.4)
            axp.plot(t, self.ana_pmean+self.ana_pstd, '--', color='dodgerblue', lw=1)
            axp.plot(t, self.ana_pmean-self.ana_pstd, '--', color='dodgerblue', lw=1)
        axp.legend(loc=0, fontsize=20)
        axp.tick_params(axis='both', labelsize=16)
        axp.set_xlabel(xlabel='t (s)', fontsize=20)
        axp.set_ylabel(ylabel='P (s)', fontsize=20)
        axp.grid(True)
        axp.set_title( 'closeorder: ' + str(self.closeorder), fontsize=20)

        axm.plot(t, self.ana_mmean, '-', color='blue', lw=3, label='model')
        if self.closeorder != 0:
            axm.fill_between(t, self.ana_mmean + self.ana_mstd,
                             self.ana_mmean - self.ana_mstd, color='dodgerblue', alpha=0.4)
            axm.plot(t, self.ana_mmean+self.ana_mstd, '--', color='dodgerblue', lw=1)
            axm.plot(t, self.ana_mmean-self.ana_mstd, '--', color='dodgerblue', lw=1)
        axm.legend(loc=0, fontsize=20)
        axm.tick_params(axis='both', labelsize=16)
        axm.set_xlabel(xlabel='t (s)', fontsize=20)
        axm.set_ylabel(ylabel='M (s)', fontsize=20)
        axm.grid(True)
        axm.set_title( 'closeorder: ' + str(self.closeorder), fontsize=20)
        if not ax:
            return figp, figm

    def plot_kmc(self, ax=None, each_sample=False):
        """ Plot kMC result only.
            Return two figure objects if parameter \'ax\' is not provided.
            """
        if ax:
            axp, axm = ax[0], ax[1]
            figp, figm = None, None
        else:
            figp, axp = plt.subplots(figsize=(8, 6))
            figm, axm = plt.subplots(figsize=(8, 6))
        
        t = self.kmc_t
        if each_sample:
            self._plot_kmcsamples([axp, axm], t, self.kmc_p, self.kmc_m)
        
        axp.plot(t, self.kmc_pmean, '-.', color='red', lw=2.5, label='kmc')
        axp.fill_between(t, self.kmc_pmean + self.kmc_pstd, self.kmc_pmean - self.kmc_pstd, color='tomato', alpha=0.5)
        axp.plot(t, self.kmc_pmean + self.kmc_pstd, '--', color='tomato', lw=0.8)
        axp.plot(t, self.kmc_pmean - self.kmc_pstd, '--', color='tomato', lw=0.8)
        axp.legend(loc=0, fontsize=20)
        axp.tick_params(axis='both', labelsize=16)
        axp.set_xlabel(xlabel='t (s)', fontsize=20)
        axp.set_ylabel(ylabel='P (s)', fontsize=20)
        axp.grid(True)
        axp.set_title( str(self.para['nsample']) + ' samples', fontsize=20)
        
        axm.plot(t, self.kmc_mmean, '-.', color='red', lw=2.5, label='kmc')
        axm.fill_between(t, self.kmc_mmean + self.kmc_mstd, self.kmc_mmean - self.kmc_mstd, color='tomato', alpha=0.5)
        axm.plot(t, self.kmc_mmean+self.kmc_mstd, '--', color='tomato', lw=0.8)
        axm.plot(t, self.kmc_mmean-self.kmc_mstd, '--', color='tomato', lw=0.8)
        axm.legend(loc=0, fontsize=20)
        axm.tick_params(axis='both', labelsize=16)
        axm.set_xlabel(xlabel='t (s)', fontsize=20)
        axm.set_ylabel(ylabel='M (s)', fontsize=20)
        axm.grid(True)
        axm.set_title( str(self.para['nsample']) + ' samples', fontsize=20)
        if not ax:
            return figp, figm

    def _plot_kmcsamples(self, ax, t, p, m):
        axp, axm = ax[0], ax[1]
        # plot each sample
        for isample in range(0, self.para['nsample']):
            axp.plot(t, p[isample], ls='-', color='gray', lw=0.6, alpha=0.6)
            axm.plot(t, m[isample], ls='-', color='gray', lw=0.6, alpha=0.6)
