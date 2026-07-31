import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from Signal import Signal
from OutputSignals import OutputSignals
from Filter import Filter
from FilterBank import FilterBank
import helpers

class Cochlea():
  '''
  Model of Cochlea
  '''
  def __init__(self, species=None, type=None, CF0=20, l_factor=3.8, length=20, xs=None, cfs=None, num_filters=None, rho=1.000, Ap=None, bp=None, Bu=None, gain_const=None, peak_magndb=None, Bpeak=None, fpeak=None, phiaccum=None, Nbeta=None, Nf=None, Qerb=None, ERBbeta=None, ERBf=None, Qn=None, Qn2=None, BWndBbeta=None, BWndBf=None, BWn2dBbeta=None, BWn2dBf=None, Sbeta=None, Sf=None, n=10, n2=3, betas=None, freqs=None):
    '''
    Initializes Cochlea. Most arguments are the same as for `FilterBank` object.

    NOTE: pressure and velocity formulas only make most sense when filter constants vary slowly

    Arguments:
      species: fast way to initialize Cochlea for certain species
      num_filters: number of filters to model cochlea with. Default is 4.
      CF0: characteristic frequency (kHz) of base of cochlea
      length: length of cochlea (mm)
      l_factor: constant factor for cochlear model (mm)
      xs: positions of filters (mm) along cochlea
      cfs: characteristic frequencies (kHz) of filters along cochlea
      rho: density of cochlear fluid (g/cm^3)
      Other arguments are defined identically to those in `Filter`/`FilterBank` \
        and the restrictions on lengths of vectors are the same as for `FilterBank`
    '''
    if species is not None:
      if num_filters is None:
        num_filters = 4
      CF0, l_factor, length = self._given_species(species)
    self.cochlea_length = length
    self.cf = (lambda x: CF0*np.exp(-x/l_factor))
    self.cf2x = (lambda cf: -l_factor*np.log(cf/CF0))

    type = 'pressure' if type is None else type

    if species is not None:
      xs = np.linspace(0, length, 4) # let user set?
      cfs = self.cf(xs)
      Ap = 0.3768 * np.exp(-0.1366 * cfs) # seems improbable
      bp = [1., 1., 1., 1.]
      Bu = 3.714 * np.exp(0.03123 * cfs) # same here
      args = {'Ap':Ap, 'bp':bp, 'Bu':Bu, 'cf':cfs}
    else:
      args = {'Ap':Ap, 'bp':bp, 'Bu':Bu, 'gain_const':gain_const, 'peak_magndb':peak_magndb,
              'Bpeak':Bpeak, 'fpeak':fpeak, 'phiaccum':phiaccum, 'Nbeta':Nbeta, 'Nf':Nf, 'Qerb':Qerb, 'ERBbeta':ERBbeta, 'ERBf':ERBf, 'Qn':Qn, 'Qn2':Qn2, 'BWndBbeta':BWndBbeta, 'BWndBbeta':BWndBf, 'BWn2dBbeta':BWn2dBbeta, 'BWn2dBf':BWn2dBf, 'Sbeta':Sbeta, 'Sf':Sf, 'n':n, 'n2':n2}
      num_filters = 1
      for v in args.values():
        if np.ndim(v) >= 1:
          num_filters = max(num_filters, len(v))

      if cfs is None:
        if xs is None:
          xs = np.linspace(0, length, num_filters)
        cfs = [self.cf(x) for x in xs]
      else:
        if xs is not None:
          raise Exception('Please provide exactly one of either a list of all locations along cochlea (xs) or a list of all characteristic frequencies (cfs)')
        xs = [-np.log(cf/CF0)*l_factor for cf in cfs]

      self.xs = xs
      args['cf'] = cfs
      args['betas'] = [betas for _ in range(num_filters)]
      args['freqs'] = [freqs for _ in range(num_filters)]
    # using filterbank here purely to automatically do all the calculations
    if type=='velocity':
      args['Bu'] = [x+1 for x in args['Bu']]
      self.filterbank = FilterBank(topology='series', type='V', **args)
    elif type=='pressure':
      self.filterbank = FilterBank(topology='series', type='P', **args)
    else:
      raise Exception('Cochlea type must be either P or V')

    apexmost_filter = self.filterbank.filters[-1]

    Ap_apex = apexmost_filter.get_consts()['Ap']
    bp_apex = apexmost_filter.get_consts()['bp']
    self.bp_apex = bp_apex
    Bu_apex = apexmost_filter.get_consts()['Bu']

    allparams = [[], [], []]
    for fil in self.filterbank.filters:
      params = fil.get_consts()
      allparams[0] += [params['Ap']]
      allparams[1] += [params['bp']]
      allparams[2] += [params['Bu']]
    self.Ap_fun = (lambda x: np.exp(np.interp(x, self.xs, np.log(np.array(allparams[0])))))
    self.bp_fun = (lambda x: np.exp(np.interp(x, self.xs, np.log(np.array(allparams[1])))))
    self.Bu_fun = (lambda x: np.exp(np.interp(x, self.xs, np.log(np.array(allparams[2])))))

    p = 1j*bp_apex - Ap_apex
    self.wavenumber = (lambda beta: (beta/l_factor) * 2 * Bu_apex * (1j*beta + Ap_apex) / ((1j*beta - p)*(1j*beta - p.conjugate())))
    self.k = self.wavenumber
    # k can also have a two variable version
    self.impedance_over_2picfx = (lambda beta: -2j * rho * beta / self.wavenumber(beta))
    self.Z_norm = self.impedance_over_2picfx
    self.impedance = (lambda beta, cf: self.impedance_over_2picfx(beta) * 2 * np.pi * cf)
    self.Z = self.impedance
    self.Y_norm = (lambda beta: 1/self.Z_norm(beta))
    self.Y = (lambda beta, cf: 1/self.Z(beta, cf))

  @classmethod
  def five_param(cls, type=None, aAp=None, bAp=None, bp=None, aBu=None, bBu=None, gain_const=None, peak_magndb=None, CF0=20, l_factor=3.8, length=20, xs=None, cfs=None, rho=1.000, betas=None, freqs=None):
    '''
    Five parameter parameterization of Cochlea from (Alkhairy 2019)
    '''
    cf = (lambda x: CF0*np.exp(-x/l_factor))
    Ap_func = (lambda x: aAp*np.exp(bAp*cf(x)))
    Bu_func = (lambda x: aBu*np.exp(bBu*cf(x)))
    if cfs is not None:
      if xs is not None:
        raise Exception('Please provide exactly one of either a list of all locations along cochlea (xs) or a list of all characteristic frequencies (cfs)')
      xs = [-np.log(cf/CF0)*l_factor for cf in cfs]
      cochlea = cls(type=type, Ap=[Ap_func(x) for x in xs], bp=bp, Bu=[Bu_func(x) for x in xs], gain_const=gain_const, peak_magndb=peak_magndb, CF0=CF0, l_factor=l_factor, length=length, cfs=cfs, rho=rho, species=None, betas=betas, freqs=freqs)
    else:
      if xs is None:
        xs = np.linspace(0, length, 4)
      cochlea = cls(type=type, Ap=[Ap_func(x) for x in xs], bp=bp, Bu=[Bu_func(x) for x in xs], gain_const=gain_const, peak_magndb=peak_magndb, CF0=CF0, l_factor=l_factor, length=length, xs=xs, rho=rho, species=None, betas=betas, freqs=freqs)
    cochlea.Ap_fun = Ap_func
    cochlea.Bu_fun = Bu_func
    return cochlea

  def filter_at_location(self, x_coord, gain_const=None, peak_magndb=None, type='P'):
    '''
    Returns filter equivalent to behavior of Cochlea at position x_coord

    If both gain_const or peak_magndb is not provided, the values are those that are provided for
    the filter nearest the base of the cochlea (i.e. if gain_const was specified for that base filter,
    gain_const will be the same for this filter. If peak_magndb was specified, peak_magndb will be
    the same for this filter).

    Arguments:
      x_coord: in mm
      type: type of Filter ('P' or 'V')
    '''
    if gain_const is None and peak_magndb is None:
      gain_const = self.filterbank.filters[0].filter.gain_const
      peak_magndb = self.filterbank.filters[0].filter.peak_magndb
    return Filter(type=type, Ap=self.Ap_fun(x_coord), bp=self.bp_fun(x_coord), Bu=self.Bu_fun(x_coord), gain_const=gain_const, peak_magndb=peak_magndb, cf=self.cf(x_coord))

  def _given_species(self, species):
    if species == 'chinchilla': # Muller et al HR 2010
      CF0 = 28.131
      l_factor = 3.6649
      length = 20
    elif species == 'human':
      CF0 = 20.823
      l_factor = 7.2382
      length = 35
    elif species == 'guinea pig' or species == 'guineapig': # Tsuji and Liberman 1997 J. Comp. Neurol. 381:188-202
      CF0 = 54.732
      l_factor = 3.3178
      length = 20
    elif species == 'mouse':
      CF0 = 71.130
      l_factor = 1.8566
      length = 5.13
    else:
      raise Exception(f'"{species}" is an unsupported species')

    return (CF0, l_factor, length)

  def plot(self, var, betas=None, cfs=None, setting='realimag', custom_title='Function plot', show=True, phase_in_rad=True):
    '''
    Plots Cochlear variables in various ways.

    Arguments:
      var: function (such as self.k, self.Z_norm)
      betas: normalized frequencies to evaluate normalized impedance at
      setting: one of 'magnphase', 'realimag', 'nichols', 'nyquist'. \
        This specifies the type of graph used to plot the normalized impedance
      custom_title: Optional title of plot. Default is 'Normalized impedance (Z_norm)'.
      show: `True` if plot is to be shown, `False` otherwise. Default is `True`.
      phase_in_rad: Show phase in radians if True or in cycles otherwise
    '''
    if cfs is None:
      if betas is None:
        betas = np.linspace(0.01, self.bp_apex*1.5, 10000)
      else:
        betas = np.array(betas) # seems like this line is unnecessary?
      fdata = var(betas)
    else:
      fdata = var(betas, cfs)
    reals = np.array([x.real for x in fdata])
    imags = np.array([x.imag for x in fdata])
    magns = abs(fdata)
    if phase_in_rad:
      phases = np.unwrap(np.angle(fdata))
    else:
      phases = np.unwrap(np.angle(fdata)) / (2*np.pi)

    if show:
      if setting == 'realimag':
        plt.plot(betas, reals)
        plt.plot(betas, imags)
        plt.axhline(y=0, color='k', ls=':')
        plt.axvline(x=self.bp_apex, color='k', ls=':')
        plt.title(custom_title)
        plt.xlabel('Normalized frequency (Hz)')
        plt.ylabel('Y_norm (???/Ω)')
        plt.show()
      if setting == 'magnphase':
        helpers.plot_2x1(betas, magns, phases, xlabel='Normalized frequency (Hz)', upper_ylabel='Magnitude(Y_norm)', lower_ylabel=f'Phase(Y_norm) ({"rad" if phase_in_rad else "cyc"})', custom_title=custom_title)
      if setting == 'nichols':
        helpers.plot_with_arrow(phases, magns, xlabel=f'Phase ({"rad" if phase_in_rad else "cyc"})', ylabel='Magnitude(Y_norm)', custom_title=custom_title)
      if setting == 'nyquist':
        helpers.plot_with_arrow(reals, imags, xlabel='Re (???/Ω)', ylabel='Im(Y_norm) (???/Ω)', custom_title=custom_title)
    return [betas, reals, imags, magns, phases]

  def signal_response_heatmap(self, signal, len_xs=20, custom_title='Envelope of Signal Response at Various Characteristic Frequencies', show=True):
    '''
    Heatmap of Cochlear response to Signal. Successively more apical parts \
      of the Cochlea are checked to see their response to the same Signal, \
      which neatly shows which parts of the Cochlea are most response to what parts of the Signal

    NOTE: errors if Cochlea has something of length 1?

    Arguments:
      signal: Signal that Cochlea is processing
      len_xs: Gives number of points along Cochlea to filter Signal at. Default is 20.
      custom_title: Optional title of plot. Default is 'Signal Heatmap'.
      show: `True` if plot is to be shown, `False` otherwise. Default is `True`.
    '''
    sigs = []
    cfs = []
    for x in np.linspace(0.025, self.cochlea_length, len_xs):
      fil = self.filter_at_location(x)
      cfs += [round(self.cf(x), 2)]
      sig = fil.solve(signal, method='tf')
      sig = sig.envelope_analytic()[0]
      sigs += [sig]

    if show:
      fig, ax = plt.subplots()
      img = ax.imshow(sigs, cmap='viridis', aspect='auto', extent=[signal.timestamps[0], signal.timestamps[-1], cfs[-1], cfs[0]])
      fig.suptitle(custom_title)
      fig.colorbar(img)
      ax.set_xlabel('Time (ms)')
      ax.set_ylabel('Characteristic Frequency (kHz)')
      plt.show()
    return sigs

def signal_response_heatmap(self, signal, filters, custom_title='Envelope of Signal Response at Various Characteristic Frequencies', show=True):
  '''
  Heatmap of Cochlear response to Signal. Successively more apical parts \
    of the Cochlea are checked to see their response to the same Signal, \
    which neatly shows which parts of the Cochlea are most response to what parts of the Signal

  Arguments:
    signal: Signal that Cochlea is processing
    len_xs: Gives number of points along Cochlea to filter Signal at. Default is 20.
    custom_title: Optional title of plot. Default is 'Signal Heatmap'.
    show: `True` if plot is to be shown, `False` otherwise. Default is `True`.
  '''
  sigs = []
  cfs = []
  for fil in filters:
    sigs += [fil.solve(signal, method='tf')]

  if show:
    plt.yscale('symlog')
    plt.pcolormesh(signal.timestamps, cfs, sigs)
    plt.show()
  return sigs
