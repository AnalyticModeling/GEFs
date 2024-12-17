import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import helpers

from FilterType import *
from Signal import Signal

class Filter:
  '''
  Class for specifying filters. Actual filter is an attribute of Filter after some dispatch on arguments
  '''
  _num_filters_made = 0
  def __init__(self, ir=None, tf=None, coeffs=None, roots=None, type=None, Ap=None, bp=None, Bu=None, gain_const=None, peak_magndb=0, Bpeak=None, fpeak=None, phiaccum=None, Nbeta=None, Nf=None, Qerb=None, ERBbeta=None, ERBf=None, Qn=None, Qn2=None, BWndBbeta=None, BWndBf=None, BWn2dBbeta=None, BWn2dBf=None, Sbeta=None, Sf=None, n=10, n2=3, betas=None, freqs=None, cf=None):
    '''
    Initializes `Filter` object from provided inputs. There are 6 ways to \
      initialize a `Filter`. See a bit later for specifications on arguments.

    - `ir`: Generates Filter given impulse response
    - `tf`: Generates Filter given transfer function
    - `coeffs`: Uses `coeffs` as coefficients of numerator/denominator of transfer function
    - `roots`: Uses `roots` as zeros/poles of transfer function
    - `Ap`, `bp`, `Bu`: Filter parameters as defined in Alkhairy (2019)
    - `Bpeak`, `fpeak`, `phiaccum`, `Nbeta`, `Nf`, `Qerb`, `ERBbeta`, `ERBf`, \
      `Qn`, `Qn2`, `BWndBbeta`, `BWndBf`, `BWn2dBbeta`, `BWn2dBf`, `Sbeta`, `Sf`: \
      Filter characteristics in Alkhairy (2019). A transfer function parameterized \
      with `Ap`, `bp`, `Bu` is created that approximates the provided set of characteristics

    Only one way of specifying a filter should be used. All 3 parameters must be provided if \
      that is used. Three characteristics should be provided if that is used, although Bpeak will \
      default to 1 if only 2 characteristics are provided and Bpeak isn't one of them.

    Either `gain_const` or `peak_magndb` should be specified. If neither is specified \
      `peak_magndb` will default to 0. `gain_const` determines the gain constant at the \
      peak of the transfer function; `peak_magndb` determines the absolute magnitude \
      (in dB) of the peak of transfer function.

    `cf` is the characteristic frequency (kHz) of the filter (see Alkhairy (2019)); calculations \
      are done by normalizing inputs to this CF. If any inputs in terms of unnormalized frequencies \
      are provided (`fpeak`, `Nf`, `ERBf`, `BWndBf`, `BWn2dBf`, `Sf`), `cf` must be provided. \
      If inputs are only in terms of normalized frequencies, `cf` is optional (but will default to 1 \
      for implementation purposes)

    `betas`/`freqs`

    phiaccum is in radians, Nbeta/Nf is in cycles.

    All transfer functions will end up as

    Arguments:
      ir: Impulse response (function from t (s) to real data (anything but often kPa), e.g. `lambda t: np.exp(-t)`)
      tf: Transfer function (function from s = j*frequency (kHz) to real data (anything but often kPa), e.g. `lambda s: 1/(1+s)`)
      coeffs: Coefficients of rational transfer function. [[1, 2], [3, 4, 5]] corresponds to the transfer function `lambda s: (1+2*s)/(3+4*s+5*s**2)`
      roots: Zeros and poles of rational transfer function. [[1], [1+2j, 1-2j]] corresponds to the transfer function `lambda s: (s-1)/((s-(1+2j))*(s-(1-2j)))`
      type: Must be either `'P'` or `'V'`. Determines whether the filter returns functions that give pressure (`'P'`) or velocity (`'P'`)
      Ap, bp, Bu: Filter parameters. See Alkhairy (2019)
      gain_const, peak_magndb: Gain constant/peak magnitude (in dB) of Parameterized filter (and by extension, filters defined from following filter characteristics)
      Bpeak, Nbeta, ERBbeta, BWndBbeta, BWn2dBbeta, Sbeta: Filter characteristics in terms of frequencies normalized with filter's characteristic frequency (Alkhairy 2019). \
        Note `Nbeta` is in cycles.
      fpeak, Nf, ERBf, BWndBf, BWn2dBf, Sf: Filter characteristics in terms of unnormalized frequencies (Alkhairy 2019). Note Nf is in cycles.
      phiaccum, Qerb, Qn, Qn2: Filter characteristics that don't rely on whether or not the filter is in defined in terms of normalized or unnormalized characteristics. \
        Note phiaccum is in radians
      n, n2: Subscripts of Qn and Qn2 if applicable
      betas, freqs: Values of beta/f to use to computationally evaluate e.g. filter characteristics. In kHz
      cf: Characteristic frequency of filter
    '''
    self.uid = self._num_filters_made
    Filter._num_filters_made += 1

    has_tf = (tf is not None)
    has_ir = (ir is not None)
    has_coeffs = (coeffs is not None)
    has_roots = (roots is not None)
    has_params = any(param is not None for param in [Ap, bp, Bu])
    # the following is incomplete
    has_chars = any(characteristic is not None for characteristic in [Bpeak, phiaccum, Nbeta, Qerb, Qn, Qn2, Sbeta])
    if sum([has_coeffs, has_roots, has_tf, has_params, has_chars]) != 1:
      raise Exception('Exactly one filter representation should be used')

    self.in_terms_of_normalized = True
    if any(arg is not None for arg in [freqs, ERBf, BWndBf, BWn2dBf, Sf, Nf]):
      if cf is None:
        raise Exception('Please provide characteristic frequency of filter')
      self.in_terms_of_normalized = False

    self.cf = 1 if cf is None else cf # ok if cf is unspecified if everything in terms of beta

    if freqs is not None: betas = [f/self.cf for f in freqs]
    if fpeak is not None: Bpeak = fpeak/self.cf
    if ERBf is not None: ERBbeta = ERBf/self.cf
    if BWndBf is not None: BWndBbeta = BWndBf/self.cf
    if BWn2dBf is not None: BWn2dBbeta = BWn2dBf/self.cf
    if Sf is not None: Sbeta = Sf * self.cf**2
    if Nf is not None: Nbeta = Nf * self.cf

    if betas is None: betas = np.geomspace(0.01, 10, 10000) # is there a more adaptive way to pick betas if it is not provided

    # eventually refactor by making the __init__ of the three classes deal with the logic (which it actually already mostly does)
    if has_tf:
      self.filter = Arbitrary(self.uid, tf=tf, betas=betas)
    elif has_ir:
      self.filter = Arbitrary(self.uid, ir=ir, betas=betas)
    elif has_coeffs:
      self.filter = Rational(self.uid, coeffs=coeffs, betas=betas)
    elif has_roots:
      self.filter = Rational(self.uid, roots=roots, betas=betas)
    elif has_params:
      self.filter = Parameterized(self.uid, type=(type if type is not None else 'P'), Ap=Ap, bp=bp, Bu=Bu, gain_const=gain_const, peak_magndb=peak_magndb, betas=betas)
    elif has_chars:
      self.filter = Parameterized(self.uid, type=(type if type is not None else 'P'), gain_const=gain_const, peak_magndb=peak_magndb, Bpeak=Bpeak, phiaccum=phiaccum, Nbeta=Nbeta, Qerb=Qerb, ERBbeta=ERBbeta, Qn=Qn, Qn2=Qn2, BWndBbeta=BWndBbeta, BWn2dBbeta=BWn2dBbeta, Sbeta=Sbeta, n=(10 if n is None else n), n2=(3 if n2 is None else n2), betas=betas)
    else:
      raise Exception('should be impossible to get here')

  @classmethod
  def multiband_params(cls, type='P', Ap=0.1, bp=1, Bu=3, gain_const=None, peak_magndb=0, betas=None, freqs=None, cf=None):
    '''
    Specifies a multiband parameterized filter of type `type` (by summing up
    transfer functions corresponding to filters filtering for each band)
    using `Ap`, `bp`, `Bu`, `gain_const`, `peak_magndb`

    `Ap`, `bp`, `Bu`, `gain_const`, `peak_magndb` must be able to be broadcasted to lists of the same length.

    See `__init__` for definitions and behavior of arguments
    '''
    if type is None: type = 'P'
    if Ap is None: Ap = 0.1
    if bp is None: bp = 1
    if Bu is None: Bu = 3

    Ap, bp, Bu, gain_const, peak_magndb = helpers.match_lengths(Ap, bp, Bu, gain_const, peak_magndb)

    if freqs is not None:
      if cf is None:
        raise Exception('Please provide peak frequency if freqs is specified')
      betas = [f/cf for f in freqs]
    elif betas is None:
      betas = np.geomspace(0.01, 3*max(bp), 10000)

    constituent_filters = []
    for i in range(len(bp)):
      constituent_filters += [Parameterized(uid=-1, type=type, Ap=Ap[i], bp=bp[i], Bu=Bu[i], gain_const=gain_const[i], peak_magndb=peak_magndb[i], betas=betas)]

    def tot_tf(s):
      return sum(fil.tf(s) for fil in constituent_filters)

    return cls(type=None, tf=tot_tf, betas=betas, freqs=freqs, cf=cf)

  @classmethod
  def multiband_chars(cls, type='P', gain_const=None, peak_magndb=0, Bpeak=None, phiaccum=None, Nbeta=None, Nf=None, Qerb=None, ERBbeta=None, ERBf=None, Qn=None, Qn2=None, BWndBbeta=None, BWndBf=None, BWn2dBbeta=None, BWn2dBf=None, Sbeta=None, Sf=None, n=10, n2=3, betas=None, freqs=None, cf=None):
    '''
    Specifying a multiband parameterized filter (by summing up
    transfer functions corresponding to filters filtering for
    each band) by using `gain_const`/`peak_magndb` and characteristics.

    Characteristics and `gain_const`/`peak_magndb` must be able to be broadcasted to lists of the same length.

    See `__init__` for definitions and behavior of arguments
    '''
    if any(arg is not None for arg in [freqs, ERBf, BWndBf, BWn2dBf, Sf, Nf]):
      if cf is None:
        raise Exception('Please provide peak frequency if freqs is specified')

    if freqs is not None:
      betas = [f/cf for f in freqs]
    elif betas is None:
      betas = np.geomspace(0.01, 10, 10000)
    if ERBf is not None: ERBbeta = [val/cf for val in ERBf]
    if BWndBf is not None: BWndBbeta = [val/cf for val in BWndBf]
    if BWn2dBf is not None: BWn2dBbeta = [val/cf for val in BWn2dBf]
    if Sf is not None: Sbeta = [val * cf**2 for val in Sf]
    if Nf is not None: Nbeta = [val * cf for val in Nf]

    type, gain_const, peak_magndb, Bpeak, phiaccum, Nbeta, Qerb, ERBbeta, Qn, Qn2, BWndBbeta, BWn2dBbeta, Sbeta, n, n2 = \
      helpers.match_lengths(type, gain_const, peak_magndb, Bpeak, phiaccum, Nbeta, Qerb, ERBbeta, Qn, Qn2, BWndBbeta, BWn2dBbeta, Sbeta, n, n2)

    constituent_filters = []
    for i in range(len(n)):
      constituent_filters += [Parameterized(uid=-1, type=type[i], gain_const=gain_const[i], peak_magndb=peak_magndb[i], Bpeak=Bpeak[i], phiaccum=phiaccum[i], Nbeta=Nbeta[i], Qerb=Qerb[i], ERBbeta=ERBbeta[i], Qn=Qn[i], Qn2=Qn2[i], BWndBbeta=BWndBbeta[i], BWn2dBbeta=BWn2dBbeta[i], Sbeta=Sbeta[i], n=n[i], n2=n2[i], betas=betas)]

    def tot_tf(s):
      return sum(fil.tf(s) for fil in constituent_filters)

    return cls(type=None, tf=tot_tf, betas=betas, freqs=freqs, cf=cf)

  def get_computed_chars(self):
    '''
    Getter for dictionary of all computed filter characteristics in terms of *normalized* frequencies. Keys are \
      {'Bpeak', 'Nbeta', 'phiaccum', 'Qerb', 'ERBbeta', 'Qn', 'n', 'BWndBbeta', 'Qn2', 'n2', 'BWn2dBbeta', 'Sbeta'}. \
      The values are approximated from the transfer function.
    '''
    return self.filter.chars

  def get_computed_unnormalized_chars(self):
    '''
    Getter for dictionary of computed filter characteristics in terms of *unnormalized* frequencies. Keys are \
      {'fpeak', 'Nf', 'phiaccum', 'Qerb', 'ERBf', 'Qn', 'n', 'BWndBf', 'Qn2', 'n2', 'BWn2dBf', 'Sf'}. \
      The values are approximated from the transfer function
    '''
    if self.in_terms_of_normalized:
      raise Exception('Filter was only initialized with inputs in terms of normalized frequencies')
    cs = self.filter.chars
    newcs = dict()
    for k in cs:
      if k == 'Bpeak': newcs['fpeak'] = cs['Bpeak']*self.cf
      elif k == 'ERBbeta': newcs['ERBf'] = cs['ERBbeta']*self.cf
      elif k == 'BWndBbeta': newcs['BWndBf'] = cs['BWndBbeta']*self.cf
      elif k == 'BWn2dBbeta': newcs['BWn2dBf'] = cs['BWn2dBbeta']*self.cf
      elif k == 'Sf': newcs['Sf'] = cs['Sbeta']/self.cf**2
      elif k == 'Nf': newcs['Nf'] = cs['Nbeta']/self.cf
    return newcs

  def get_orig_chars(self):
    '''
    Getter for dictionary of originally provided filter characteristics. These may differ slightly from \
      the outputs of `get_computed_chars` and `get_computed_unnormalized_chars`. For instance, if Filter is \
      initialized with Filter(Bpeak=1, Nbeta=19.1, Qerb=25.9), then this function will return \
      {'Bpeak':1, 'Nbeta':19.1, 'Qerb':25.9}. The resulting Filter that best fits these characteristics may have \
      characteristics slightly off due to the use of numerical approximations.
    '''
    if isinstance(self.filter, Parameterized):
      return self.filter.orig_chars
    raise Exception(f'Original characteristics undefined')

  def get_params(self):
    '''
    Getter for dictionary of originally provided filter parameters. For instance, if Filter is \
      initialized with Filter(Ap=0.05, bp=1, Bu=3), then this function will return \
      {'Ap':0.05, 'bp':1, 'Bu':3}
    '''
    if isinstance(self.filter, Parameterized):
      return self.filter.params
    raise Exception(f'Parameters undefined')

  def solve(self, input, method=None, fs=None):
    '''
    Finds output signal given input signal. There are five allowed values of `method`:

    - None: If `method` is unspecified, a default method will be picked (currently uses 'tf' \
        everywhere except for `ir` when the Filter was originally instantiated using 'ir')
    - 'tf': Multiplies signal in frequency domain by corresponding values of transfer function and FFTs to get output signal
    - 'ir': Convolves signal in time domain by impulse response
    - 'ode': If Filter is Rational OR Parameterized with `Bu` an integer, the transfer function will be a rational function, \
        in which case there will be a corresponding ODE that can be solved numerically to get the output signal. If `Bu` is not \
        and integer, `Bu` will be rounded to an integer and there will be a warning.
    - 'fde': If Filter is Parameterized, regardless of `Bu`, the transfer function will correspond to a fractional differential \
        equation, which can be solved numerically with a Riemann-Liouville integral.

    If method is 'ode', a higher sampling rate of input data may be required.

    If `input` is not a `Signal` object, the sampling frequency of the input should be specified
    with `fs` (although it will default to 1 if undefined).

    If `input` is not a `Signal` object, `input` is assumed to be in terms of real time, not scaled time.

    If `input` is/isn't a `Signal` object, the output will/won't be a `Signal` object.

    Attributes:
      input: input signal
      method: method used to calculate output signal. Must be one of (None, 'tf', 'ir', 'ode', 'fde').
      fs: Sampling frequency of input (if input is not a Signal object)
    '''
    if isinstance(input, Signal):
      if input.mode in ['t', 'f', 'w']:
        sig = Signal(mode='ttilde', data=input.mode_t, fs=input.fs/(2*np.pi*self.cf)) # <- check
      else:
        sig = input
    else:
      if fs is None:
        raise Exception('Must provide fs')
      sig = Signal(mode='ttilde', data=input, fs=fs/(2*np.pi*self.cf))

    if method is None:
      if isinstance(self.filter, Arbitrary) and not self.filter.init_with_tf:
        method = 'ir'
      else:
        method = 'tf'

    if method == 'tf':
      output = self.filter.solve_fourier_tf(len(sig), sig['b'], sig.freqstamps)
    elif method == 'ir':
      # if isinstance(self.filter, Parameterized):
      #   warnings.warn('Parametrized approximate impulse response may produce incorrect values for very small t')
      output = self.filter.solve_convolve_ir(sig['ttilde'], sig.timestamps)
    elif method == 'ode':
      if isinstance(self.filter, Rational):
        output = self.filter.solve_diffeq(sig['ttilde'], fs=sig.fs)
      elif isinstance(self.filter, Parameterized):
        output = self.filter.solve_as_rational_tf(sig['ttilde'], fs=sig.fs)
      else:
        raise Exception(f'Numerical integration of transfer function not available for non-rational or non-parameterized filters')
    elif method == 'fde':
      if not isinstance(self.filter, Parameterized):
        raise Exception(f'Fractional integration of transfer function not available for non-parameterized filters')
      output = self.filter.solve_fractional_diffeq(sig['ttilde'], fs=sig.fs)
    else:
      raise Exception('Invalid solution method')

    if isinstance(input, Signal) and input.mode in ['b', 'beta', 'ttilde']:
      return Signal(mode='ttilde', data=output, fs=sig.fs)
    else:
      return Signal(mode='t', data=output, fs=sig.fs*2*np.pi*self.cf)

  def __call__(self, input, method='default', fs=1):
    return self.solve(input, method=method, fs=fs)

  def bode_plot(self, freqs=None, peak_magndb=0, custom_title='Bode plot', show=True):
    '''
    Generate Bode plot of filter. Returns [x-axis (frequency) data, magnitudes (dB), phases (cycles)].

    Inaccurate results may occur if `freqs` is not wide enough.

    Note transfer functions for Filters are internally stored as functions of normalized frequency, so `freqs` should also be normalized frequencies.

    Attributes:
      freqs: Frequencies (kHz) where transfer function is sampled (after multiplication by 1j)
      peak_magndb: Peak magnitude that should be displayed on Bode plot (defaults to 0)
      custom_title: Optional title of plot. Default is 'Bode plot'.
      show: `True` if plot is to be shown, `False` otherwise. Default is `True`.
    '''
    if freqs is None:
      xaxis = np.geomspace(0.1, 3*self.filter.chars['Bpeak'], 10000)
    else:
      xaxis = np.array(freqs)

    response = self.filter.tf(1j*xaxis)
    magns = helpers.mag2db(abs(response))
    if peak_magndb is not None:
      magns += peak_magndb-max(magns)
    phases = np.unwrap(np.angle(response)) / (2 * np.pi)

    if show:
      fig, (ax1, ax2) = plt.subplots(2, 1)
      fig.suptitle(custom_title)

      ax1.semilogx(xaxis, magns) # magn in db
      ax1.xaxis.set_major_locator(locator=matplotlib.ticker.LogLocator(subs=(0.5, 1, 2)))
      ax1.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
      ax1.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
      ax1.set_ylabel('Magnitude (dB)')

      ax2.semilogx(xaxis, phases) # phase in cycles
      ax2.xaxis.set_major_locator(locator=matplotlib.ticker.LogLocator(subs=(0.5, 1, 2)))
      ax2.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
      ax2.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
      ax2.set_ylabel('Phase (cycles)')
      ax2.set_xlabel('Normalized frequency')

      plt.show()
    return [xaxis, magns, phases]

  def frequency_real_imag_plot(self, freqs=None, custom_title='Frequency response plot', show=True):
    '''
    Similar to `bode_plot` except graphs real and imaginary parts of transfer function. Return [x-axis (frequency) data, real parts (kPa), imaginary parts (kPa)].

    Inaccurate results may occur if `freqs` is not wide enough.

    Note transfer functions for Filters are internally stored as functions of normalized frequency, so `freqs` should also be normalized frequencies.

    Attributes:
      freqs: Frequencies (kHz) where transfer function is sampled (after multiplication by 1j)
      custom_title: Optional title of plot. Default is 'Frequency response plot'.
      show: `True` if plot is to be shown, `False` otherwise. Default is `True`.
    '''
    if freqs is None:
      xaxis = np.geomspace(0.1, 1.5*self.filter.chars['Bpeak'], 10000)
    else:
      xaxis = np.array(freqs)

    response = self.filter.tf(1j*xaxis)
    reals = [z.real for z in response]
    imags = [z.imag for z in response]

    if show:
      fig, (ax1, ax2) = plt.subplots(2, 1)
      fig.suptitle(custom_title)

      ax1.semilogx(xaxis, reals)
      ax1.xaxis.set_major_locator(locator=matplotlib.ticker.LogLocator(subs=(1, 2, 4, 6, 8)))
      ax1.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
      ax1.set_ylabel('Re(z)')

      ax2.semilogx(xaxis, imags)
      ax2.xaxis.set_major_locator(locator=matplotlib.ticker.LogLocator(subs=(1, 2, 4, 6, 8)))
      ax2.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
      ax2.set_ylabel('Im(z)')
      ax2.set_xlabel('Normalized frequency')

      plt.show()
    return [xaxis, reals, imags]

  def nichols_plot(self, freqs=None, custom_title='Nichols plot', show=True):
    '''
    Generate Nichols plot of Filter. Return [x-axis (frequency) data, phases (cycles), magnitudes (dB)].

    NOTE: phases and magnitudes are swapped in output order relative to Bode plot

    Note transfer functions for Filters are internally stored as functions of normalized frequency, so `freqs` should also be normalized frequencies.

    Attributes:
      freqs: Frequencies (kHz) where transfer function is sampled (after multiplication by 1j)
      custom_title: Optional title of plot. Default is 'Nichols plot'.
      show: `True` if plot is to be shown, `False` otherwise. Default is `True`.
    '''
    if freqs is None:
      xaxis = np.geomspace(0.1, 2*self.filter.chars['Bpeak'], 10000)
    else:
      xaxis = np.array(freqs)
    normalized_magn_tf = (lambda s: self.filter.tf(s)/self.filter.tf(1j*self.filter.chars['Bpeak']))
    response = normalized_magn_tf(1j*xaxis)

    phases = np.unwrap(np.angle(response))
    magns = helpers.mag2db(abs(response))

    if show:
      plt.plot(phases, magns)
      plt.title(custom_title)
      plt.xlabel('Phase (rad)')
      plt.ylabel('Magnitude (dB)')
      # pick better arrows?
      # arrow_length = 2
      # dx = phases[-1]-phases[-2]
      # dy = magns[-1]-magns[-2]
      # theta = np.arctan2(dy, dx)
      # plt.arrow(phases[-2], magns[-2], np.cos(theta)*arrow_length, np.sin(theta)*arrow_length, head_width=0.3, head_length=3, length_includes_head=True)
      plt.show()
    return [xaxis, phases, magns]

  def nyquist_plot(self, freqs=None, custom_title='Nyquist plot', show=True):
    '''
    Generate Nyquist plot of Filter. Return [x-axis (frequency) data, real parts (kPa), imaginary parts (kPa)].

    Note transfer functions for Filters are internally stored as functions of normalized frequency, so `freqs` should also be normalized frequencies.

    Attributes:
      freqs: Frequencies (kHz) where transfer function is sampled (after multiplication by 1j)
      custom_title: Optional title of plot. Default is 'Nyquist plot'.
      show: `True` if plot is to be shown, `False` otherwise. Default is `True`.
    '''
    if freqs is None:
      xaxis = np.geomspace(0.1, 2*self.filter.chars['Bpeak'], 10000)
    else:
      xaxis = np.array(freqs)
    normalized_magn_tf = (lambda s: self.filter.tf(s)/self.filter.tf(1j*self.filter.chars['Bpeak']))
    response = normalized_magn_tf(1j*xaxis)
    reals = np.real(response)
    imags = np.imag(response)
    if show:
      plt.plot(reals, imags)
      plt.title(custom_title)
      plt.xlabel('Re(z)')
      plt.ylabel('Im(z)')
      # play around with setting
      # arrow_length = 0.05
      # dx = reals[-1]-imags[-2]
      # dy = reals[-1]-imags[-2]
      # theta = np.arctan2(dy, dx)
      # plt.arrow(reals[-2], imags[-2], np.cos(theta)*arrow_length, np.sin(theta)*arrow_length, head_width=0.01, head_length=0.01, length_includes_head=True)
      plt.show()
    return [xaxis, reals, imags]

  def impulse_response_plot(self, times=None, custom_title='Impulse response', show=True):
    '''
    Generate plot of impulse reponse of Filter. Output is [x-axis (time) data, impulse response].

    Attributes:
      times: Timestamps (s) where impulse response is evaluated. Defaults to `np.linspace(0, 100, 10000)`.
      custom_title: Optional title of plot. Default is 'Nyquist plot'.
      show: `True` if plot is to be shown, `False` otherwise. Default is `True`.
    '''
    # find optimal num_samples and t automatically?
    # slow if tf is used to define; is there any way to be more efficient
    if times is None:
      xaxis = np.linspace(0, 100, 10000)
    else:
      xaxis = np.array(times)
    response = self.filter.ir(xaxis)

    if show:
      plt.plot(xaxis, response)
      plt.title(custom_title)
      plt.xlabel('Time (s)')
      plt.show()
    return [xaxis, response]

  def pole_zero_plot(self, custom_title=None, show=True):
    '''
    Generate pole-zero plot. Output is [list of zeros, list of poles], which is the same order as inputted into Filter to initalize it.

    In graph, zeros are orange o's and poles are blue x's.

    Only makes sense for Rational and Parameterized filters. Parameterized filters have pole-zero plots of their base filters (so ignoring the `Bu` exponent)

    Attributes:
      custom_title: Optional title of plot. Default is 'Nyquist plot'.
      show: `True` if plot is to be shown, `False` otherwise. Default is `True`.
    '''
    if isinstance(self.filter, Rational):
      zeros, poles = self.filter.roots
    elif isinstance(self.filter, Parameterized):
      p = -self.filter.params['Ap'] + 1j*self.filter.params['bp']
      poles = [p, p.conjugate()]
      if self.filter.type == 'V':
        zeros = [-self.filter.params['Ap']]
      else:
        zeros = []
    else:
      raise Exception(f'Pole-zero plot unavailable for filters with arbitrary transfer functions')
    if show:
      fig, ax = plt.subplots()
      if custom_title is None:
        if isinstance(self.filter, Parameterized) and abs(self.filter.params['Bu']-round(self.filter.params['Bu']))>1e-10:
          fig.suptitle('Pole-zero plot (of base filter)')
        else:
          fig.suptitle('Pole-zero plot')
      else:
        fig.suptitle(custom_title)
      ax.axhline(y=0, color='k', ls=':')
      ax.axvline(x=0, color='k', ls=':')
      ax.scatter([z.real for z in zeros], [z.imag for z in zeros], marker='o', facecolors='none', edgecolors='tab:orange')
      ax.scatter([z.real for z in poles], [z.imag for z in poles], marker='x', edgecolors='tab:blue')
      ax.set_xlabel('Re(z)')
      ax.set_ylabel('Im(z)')
      plt.axis('equal')
      plt.show()
    return [zeros, poles]

  def Qn_plot(self, max_n=20, custom_title='Qn plot', show=True):
    '''
    Generate plot of n vs Qn. n (dB) is an indication of how wide a window to consider for Q factor.

    Returns [values of n, values of Qn]

    Attributes:
      max_n: Qn will be calculated for values of n between 1 and max_n inclusive
      custom_title: Optional title of plot. Default is 'Nyquist plot'.
      show: `True` if plot is to be shown, `False` otherwise. Default is `True`.
    '''
    if not isinstance(self.filter, Parameterized):
      raise Exception('Qn plot not guaranteed')
    params = self.filter.params
    xaxis = range(1, max_n+1)
    Qns = [params['bp'] / 2 / params['Ap'] / (10 ** (n / 10 / params['Bu']) -1) ** 0.5 for n in xaxis]
    if show:
      plt.plot(xaxis, Qns)
      plt.title(custom_title)
      plt.xlabel('n (dB)')
      plt.ylabel('Q$_n$')
      plt.show()
    return [xaxis, Qns]

  def characteristic_error(self, custom_title='Estimated vs Desired Bar Chart', show=True):
    '''
    Graphs bar chart of error rate for original characteristics used to \
      initialize filter if they exist. For instance, if filter was specified \
      from Bpeak, Nbeta, and phiaccum, those values are used to calculate \
      parameters, which are then used to generate the transfer function for \
      the filter. The approximate values of Bpeak, Nbeta, and phiaccum that \
      can be calculated from that transfer function are then compared to the \
      originally specified values to compute the errors.

    Graph is invalid if filter was not specified from original characteristics.

    Returns a dictionary with keys that are the names of the original characteristics \
      and values that are the corresponding errors.

    Attributes:
      custom_title: Optional title of plot. Default is 'Estimated vs Desired Bar Chart'.
      show: `True` if plot is to be shown, `False` otherwise. Default is `True`.
    '''
    if not isinstance(self.filter, Parameterized):
      raise Exception()
    orig = self.filter.orig_chars
    if orig is None: raise Exception() # should be unneeded
    potential_errorable_chars = ['BWndBbeta', 'BWn2dBbeta', 'ERBbeta', 'Nbeta', 'Sbeta', 'Bpeak', 'phiaccum', 'Qn', 'Qn2', 'Qerb']
    errorable_chars = []
    errors = []
    for char in potential_errorable_chars:
      if char in orig:
        val = orig[char]
        errors += [self.filter.chars[char]/val-1]
        name = char
        if not self.in_terms_of_normalized:
          if char in ['BWndBbeta', 'BWn2dBbeta', 'ERBbeta', 'Sbeta', 'Nbeta']:
            name = char[:-4]+'f'
          if char == 'Bpeak':
            name = 'fpeak'
        errorable_chars += [name]

    mini = min(errors)
    maxi = max(errors)

    if show:
      fig, ax = plt.subplots()
      ax.bar(errorable_chars, errors)
      ax.set_ylim([min(mini*1.1, -0.01), max(maxi*1.1, 0.01)])
      fig.suptitle(custom_title)
      plt.show()
    return dict(zip(errorable_chars, errors))