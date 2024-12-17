import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from Signal import Signal
from OutputSignals import OutputSignals
from Filter import Filter
from FilterBank import FilterBank
import helpers

class Cochlea(FilterBank):
  def __init__(self, types=None, species=None, CF0=20, length=20, l_factor=3.8, xs=None, cfs=None, rho=1.000, Ap=None, bp=None, Bu=None, gain_const=None, peak_magndb=None, Bpeak=None, fpeak=None, phiaccum=None, Nbeta=None, Nf=None, Qerb=None, ERBbeta=None, ERBf=None, Qn=None, Qn2=None, BWndBbeta=None, BWndBf=None, BWn2dBbeta=None, BWn2dBf=None, Sbeta=None, Sf=None, n=10, n2=3, betas=None, freqs=None):
    '''
    Cochlea. All relevant properties can be calculated from beta

    Attributes:
      CF0: characteristic frequency of base of cochlea (in kHz)
      length: length of cochlea (in mm)
      l_factor: constant factor in cochlea (lowercase l)
      xs: positions of filters along cochlea
      cfs: characteristic frequencies of filters along cochlea
      rho: density of cochlear fluid
      Other arguments are defined identically to those in Filter
    '''
    # if species is not None:
    #   if species == 'chinchilla':
    #     CF0 = 28.131
    #     l_factor = 3.6649
    #     length = 35
    #   elif species == 'human':
    #     CF0 = 20.823
    #     l_factor = 7.2382
    #     length = 20
    #   # dramatic different filter at base?
    #   xs = np.linspace(0, length, 5)[1:]
    #   Ap = [0.3768 * np.exp(-0.1366 * CF0 * np.exp(-x/l_factor)) for x in xs] # this seems improbable
    #   bp = [1., 1., 1., 1.]
    #   Bu = [3.714 * np.exp(0.03123 * CF0 * np.exp(-x/l_factor)) for x in xs] # same here
    #   # print(Ap, bp, Bu)
    # else:
    #   l_factor = 3.8

    if species is not None:
      CF0, l_factor, length = self._given_species(species)

    self.cochlea_length = length

    args = {'Ap':Ap, 'bp':bp, 'Bu':Bu, 'gain_const':gain_const, 'peak_magndb':peak_magndb,
            'Bpeak':Bpeak, 'fpeak':fpeak, 'phiaccum':phiaccum, 'Nbeta':Nbeta, 'Nf':Nf, 'Qerb':Qerb, 'ERBbeta':ERBbeta, 'ERBf':ERBf, 'Qn':Qn, 'Qn2':Qn2, 'BWndBbeta':BWndBbeta, 'BWndBbeta':BWndBf, 'BWn2dBbeta':BWn2dBbeta, 'BWn2dBf':BWn2dBf, 'Sbeta':Sbeta, 'Sf':Sf, 'n':n, 'n2':n2, 'betas':betas, 'freqs':freqs}
    num_filters = 1
    for v in args.values():
      if np.ndim(v) >= 1:
        num_filters = max(num_filters, len(v))

    self.cf = (lambda x: CF0*np.exp(-x/l_factor))
    if cfs is None:
      if xs is None:
        # raise Exception('Either a list of all locations along cochlea or a list of all characteristic frequencies must be provided') # should this be necessary
        xs = np.linspace(0, length, num_filters)
      cfs = [self.cf(x) for x in xs]
    else:
      if xs is not None:
        raise Exception('Please provide either only a list of all locations along cochlea or a list of all characteristic frequencies')
      xs = [-np.log(c/CF0)*l_factor for c in cfs]
    self.xs = xs
    types = 'P' if types is None else types

    args['cf'] = cfs
    super().__init__(topology='parallel', type=types, **args)

    apexmost_filter = self.filters[-1]

    Ap_apex = apexmost_filter.get_params()['Ap']
    bp_apex = apexmost_filter.get_params()['bp']
    self.bp_apex = bp_apex
    Bu_apex = apexmost_filter.get_params()['Bu']
    # print('apex', Ap_apex, bp_apex, Bu_apex)

    self.Ap_fun = (lambda x: np.exp(np.interp(x, self.xs, np.log(np.array(Ap)))))
    self.bp_fun = (lambda x: np.exp(np.interp(x, self.xs, np.log(np.array(bp)))))
    self.Bu_fun = (lambda x: np.exp(np.interp(x, self.xs, np.log(np.array(Bu)))))

    # beta = w/2pi/CF(x)

    p = 1j*bp_apex - Ap_apex
    # k and Z both normalized to not depend on l
    self.wavenumber = (lambda beta: (beta/l_factor) * 2 * Bu_apex * (1j*beta + Ap_apex) / ((1j*beta - p)*(1j*beta - p.conjugate())))
    self.k = self.wavenumber
    self.impedance_over_2picfx = (lambda beta: -2j * rho * beta / self.wavenumber(beta))
    self.Z_norm = self.impedance_over_2picfx
    self.impedance = (lambda beta, x: self.impedance_over_2picfx(beta) * 2 * np.pi * self.cf(x))
    self.Z = self.impedance

  @classmethod
  def five_param(cls, types=None, aAp=None, bAp=None, bp=None, aBu=None, bBu=None, gain_const=None, peak_magndb=None, CF0=20, length=20, xs=None, rho=1.000):
    '''
    Five parameter parameterization of Cochlea from (Alkhairy 2019)
    '''
    if xs is None:
      xs = np.linspace(0, length, 4)
    cf = (lambda x: CF0*np.exp(-x/3.8))
    # print('cfx', cf(0))
    # print(bAp, np.exp(bAp*cf(0)))
    Ap_func = (lambda x: aAp*np.exp(bAp*cf(x)))
    Bu_func = (lambda x: aBu*np.exp(bBu*cf(x)))
    cochlea = cls(types=types, Ap=[Ap_func(x) for x in xs], bp=bp, Bu=[Bu_func(x) for x in xs], gain_const=gain_const, peak_magndb=peak_magndb, CF0=CF0, length=length, xs=xs, rho=rho, species=None)
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
      gain_const = self.filters[0].filter.gain_const
      peak_magndb = self.filters[0].filter.peak_magndb
    return Filter(type=type, Ap=self.Ap_fun(x_coord), bp=self.bp_fun(x_coord), Bu=self.Bu_fun(x_coord), gain_const=gain_const, peak_magndb=peak_magndb, cf=self.cf(x_coord))

  def _given_species(self, species):
    if species is not None:
      if species == 'chinchilla':
        CF0 = 28.131
        l_factor = 3.6649
        length = 35
      elif species == 'human':
        CF0 = 20.823
        l_factor = 7.2382
        length = 20
      # elif species == 'guinea pig' or species == 'guineapig':
      #   CF
      # dramatic different filter at base?
      xs = np.linspace(0, length, 5)[1:]
      Ap = [0.3768 * np.exp(-0.1366 * CF0 * np.exp(-x/l_factor)) for x in xs] # this seems improbable
      bp = [1., 1., 1., 1.]
      Bu = [3.714 * np.exp(0.03123 * CF0 * np.exp(-x/l_factor)) for x in xs] # same here
      # print(Ap, bp, Bu)
    else:
      l_factor = 3.8

    return CF0, l_factor, length

  def plot_wavenumber(self, betas=None, setting='realimag', custom_title='Wavenumber (k)', show=True, phase_in_rad=True):
    '''
    Plot wavenumber function of Cochlea in various ways.

    Arguments:
      betas: normalized frequencies to evaluate wavenumber at
      setting: one of 'magnphase', 'realimag', 'nichols', 'nyquist'. \
      This specifies the type of graph used to plot the wavenumber
      custom_title:
      show:
      phase_in_rad: Show phase in radians if True or in cycles otherwise
    '''
    # default of rad/cyc seems to be different between magnphase and nichols?

    # how to use params2chars
    if setting not in ['magnphase', 'realimag', 'nichols', 'nyquist']:
      raise Exception("Setting must be one of 'magnphase', 'realimag', 'nichols', 'nyquist'")
    if custom_title is None:
      custom_title = 'Wavenumber (k)'
    if betas is None:
      betas = np.linspace(0.01, self.bp_apex*1.5, 10000)
    kdata = self.k(betas)
    reals = np.real(kdata)
    imags = np.imag(kdata)
    magns = abs(kdata)
    if phase_in_rad:
      phases = np.unwrap(np.angle(kdata))
    else:
      phases = np.unwrap(np.angle(kdata)) / (2*np.pi)

    if show:
      if setting == 'realimag':
        plt.plot(betas, reals)
        plt.plot(betas, imags)
        plt.axhline(y=0, color='k', ls=':')
        plt.axvline(x=self.bp_apex, color='k', ls=':')
        plt.title(custom_title)
        plt.xlabel('Normalized frequency (Hz)')
        plt.ylabel('k (1/mm)')
        plt.show()
      if setting == 'magnphase':
        helpers.plot_2x1(betas, magns, phases, xlabel='Normalized frequency (Hz)', upper_ylabel='Magnitude(k)', lower_ylabel=f'Phase(k) ({"rad" if phase_in_rad else "cyc"})', custom_title=custom_title)
      if setting == 'nichols':
        helpers.plot_with_arrow(phases, magns, xlabel=f'Phase(k) ({"rad" if phase_in_rad else "cyc"})', ylabel='Magnitude(k)', custom_title=custom_title)
      if setting == 'nyquist':
        helpers.plot_with_arrow(reals, imags, xlabel='Re(k) (1/mm)', ylabel='Im(k) (1/mm)', custom_title=custom_title)
    return [betas, reals, imags, magns, phases]

  def plot_impedance(self, betas=None, setting='realimag', custom_title='Normalized impedance (Z_norm)', show=True, phase_in_rad=True):
    '''
    Almost identical to plot_wavenumber except all graphs are generated with Z_norm not wavenumber
    '''
    # plot Z and normalized Z
    if betas is None:
      betas = np.linspace(0.01, self.bp_apex*1.5, 10000)
    Zdata = self.Z_norm(betas)
    reals = np.real(Zdata)
    imags = np.imag(Zdata)
    magns = abs(Zdata)
    if phase_in_rad:
      phases = np.unwrap(np.angle(Zdata))
    else:
      phases = np.unwrap(np.angle(Zdata)) / (2*np.pi)

    if show:
      if setting == 'realimag':
        plt.plot(betas, reals)
        plt.plot(betas, imags)
        plt.axhline(y=0, color='k', ls=':')
        plt.axvline(x=self.bp_apex, color='k', ls=':')
        plt.title(custom_title)
        plt.xlabel('Normalized frequency (Hz)')
        plt.ylabel('Z_norm (Ω/???)')
        plt.show()
      if setting == 'magnphase':
        helpers.plot_2x1(betas, magns, phases, xlabel='Normalized frequency (Hz)', upper_ylabel='Magnitude(Z_norm)', lower_ylabel=f'Phase(Z_norm) ({"rad" if phase_in_rad else "cyc"})', custom_title=custom_title)
      if setting == 'nichols':
        helpers.plot_with_arrow(phases, magns, xlabel=f'Phase(Z_norm) ({"rad" if phase_in_rad else "cyc"})', ylabel='Magnitude(Z_norm)', custom_title=custom_title)
      if setting == 'nyquist':
        helpers.plot_with_arrow(reals, imags, xlabel='Re(Z_norm) (Ω/???)', ylabel='Im(Z_norm) (Ω/???)', custom_title=custom_title)
    return [betas, reals, imags, magns, phases]

  def signal_response_heatmap(self, signal):
    '''
    Heatmap of Cochlear response to Signal. Successively more apical part \
    of the Cochlea are checked to see their response to the same Signal
    '''
    def solvehere(data, fs, Ap, bp, Bu):
      p = -Ap + 1j*bp
      tf = lambda s: ((s - p) * (s - p.conjugate()))**(-Bu)
      X = sp.fft.rfft(data)
      freqs = sp.fft.rfftfreq(len(data), 1/fs)
      # freqs = np.linspace(0, 1, len(X)) * fs/2 #/ self.cf(x)
      # freqs /= max(abs(freqs))
      # plt.plot(freqs, abs(tf(1j*freqs)))
      # plt.title('bodelike')
      # plt.show()
      tf_vals = tf(1j*freqs)
      tf_vals /= max(abs(tf_vals))
      Y = tf_vals*X
      y = sp.fft.irfft(Y, len(data))
      return y

    def paramhere(x):
      return {'Ap':0.3768*np.exp(-0.1366*self.cf(x)), 'bp':1, 'Bu':3.714*np.exp(-0.03123*self.cf(x))}

    sigs = []
    cfs = []
    # for x in np.linspace(0.025, self.cochlea_length, 30):
    #   fil = self.filter_at_location(x)
    #   # print(fil.get_params()['Ap'])
    #   # print(fil.get_computed_chars()['Bpeak'])
    #   cfs += [round(fil.cf, 2)]
    #   # print(fil.cf)
    #   sig = fil.solve(signal, method='tf')
    #   # sig.plot()
    #   # sigs += [[abs(v)/maxsig for v in sig.mode_t]]
    #   sig = sig.envelope_analytic()[0]
    #   maxsig = max(sig)
    #   sigs += [[v/maxsig for v in sig]]
    # print(sigs)
    signal.plot()
    # print(self.cochlea_length)
    for x in np.linspace(0.025, self.cochlea_length, 20):
    # for x in 
      fil = self.filter_at_location(x)
      # params = fil.get_params()
      cfs += [round(self.cf(x), 2)]
      # print('fs', signal.fs)
      # sig = fil.solve(signal, method='tf')
      # print(self.cf(x), fil.cf)
      # sig = Signal(data=solvehere(signal['t'], signal.fs/fil.cf, **params))
      sig = fil.solve(signal, method='tf')
      # sig.plot(custom_title=str(self.cf(x)))
      sig = sig.envelope_analytic()[0]
      sigs += [sig]
    # print(sigs)


    fig, ax = plt.subplots()
    img = ax.imshow(sigs, cmap='viridis', aspect='auto')
    # fig.suptitle(custom_title)
    # fig.colorbar(img)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Frequency (1/s)')
    plt.yticks(range(len(cfs)), cfs)
    plt.show()


  # def scale_signal(self, filter, signal, unscale=False):
  #   '''
  #   backwards?
  #   '''
  #   if not isinstance(signal, Signal):
  #     raise Exception()
  #   scale_factor = 2*np.pi * self.cf(self.xs[filter.uid])
  #   if unscale:
  #     return Signal(signal['t'], fs=signal.fs*scale_factor)
  #   else:
  #     return Signal(signal['t'], fs=signal.fs/scale_factor)

if __name__ == "__main__":
  # print(np.linspace(0, 1, 2, endpoint=True))
  # print(np.interp(0, [1], [1]))
  # c = Cochlea(Ap=[0.3768*np.exp(0.01*i) for i in range(2)], bp=[0.5, 2], Bu=[3.714*np.exp(0.01*i) for i in range(2)])
  # c = Cochlea(Ap=[0.3768*np.exp(0.01*i) for i in range(4)], bp=[0.5, 1, 1.5, 2], Bu=[3.714*np.exp(0.01*i) for i in range(4)], length=3.8)
  # c = Cochlea(species='chinchilla')
  # c.bode_plot()
  # c.bode_plot(peak_magn=None)
  # f = c.filters[3]
  # f.bode_plot(peak_magn=None)
  # print(c.bode_plot())
  # c.plot_wavenumber(setting='realimag')
  # c.plot_wavenumber(setting='magnphase')
  # c.plot_wavenumber(setting='nichols')
  # c.plot_wavenumber(setting='nyquist')

  # c.plot_impedance(setting='realimag')
  # c.plot_impedance(setting='magnphase')
  # c.plot_impedance(setting='nichols')
  # c.plot_impedance(setting='nyquist')

  # c = Cochlea.five_param(types='V', aAp=0.3768, bAp=-0.1366, bp=[0.2, 0.5, 2, 5, 15], aBu=3.714, bBu=0.03123, xs=[i for i in range(5)])
  c = Cochlea.five_param(types='V', aAp=0.3768, bAp=-0.1366, bp=[1, 1, 1, 1, 1], aBu=3.714, bBu=0.03123, xs=[i for i in range(5)])
  # func = lambda t: np.exp(-1/15*(t-20)**2) * np.cos(t)
  # sig = Signal(data=[func(t) for t in range(4000)], fs=40)
  pairs = [(1.5, 200), (8, 400), (1.5, 700), (0.3, 400)]
  def tones(t):
    ans = 0
    for i in range(4):
      fi, ti = pairs[i]
      ans += np.exp(-((t-ti)/50)**2) * np.sin(2*np.pi*fi*t)
    return ans
  sig = Signal(mode='t', data=[tones(t/100) for t in range(100000)], fs=100)
  # sig.plot()
  # fil = c.filter_at_location(3)
  # print(fil.get_params())
  # fil.bode_plot(freqs=np.linspace(2, 10, 10000))

  # newfil = Filter(Ap=0.01, bp=1, Bu=3, cf=.6)
  # newfil.bode_plot()
  # newfil.solve(sig, method='ode').plot()
  c.signal_response_heatmap(sig)



  # example = 'red_truth'
  # samplerate, data = sp.io.wavfile.read(f'test_signals/{example}.wav')
  # print('example sr:', samplerate)
  # resample_factor = data.shape[0]//10_000
  # samplerate /= resample_factor
  # if len(data.shape) == 1:
  #   mono = data
  # else:
  #   mono = [sum(d) for d in data]
  # mono = [mono[i] for i in range(0, len(mono), resample_factor)]
  # mono = mono[:1000]
  # s = Signal(mode='t', data=mono, fs=samplerate)
  # s.plot()

  # s = Signal.linear_chirp(f_init=1, f_final=3, fs=1000, n=2000)
  # osigs = c.process_signal(s)
  # # print(osigs)
  # print(len(osigs))
  # for s in osigs.outsignals[3:]:
  #   plt.plot(s.timestamps, s['t'])
  # # s1 = osigs.outsignals[0]
  # plt.show()

  # f.solve(s).plot()

  ###

  # osigs = c.process_signal(s)
  # for sig in osigs.signals:
  #   sig /= max(sig['t'])
  #   # yaxis = [abs(v) for v in sig['f']]
  #   # plt.plot(range(len(sig['f'])), yaxis)
  #   plt.plot(sig.timestamps, sig['t'])
  # # osigs.correlogram()

  # # for fil in c.filters:
  # #   if not isinstance(fil, Filter): raise Exception()
  # #   sig = fil.solve(s)
  # #   sig /= max(sig['t'])
  # #   plt.plot(sig.timestamps, sig['t'])
  # plt.show()

  # c.plot_wavenumber()
  # c.plot_impedance()