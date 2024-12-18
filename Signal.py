import numpy as np
import scipy as sp
import scipy.fft
import matplotlib.pyplot as plt

def tolerate(x, eps=1e-10):
  return abs(x-round(x)) <= eps

class Signal:
  '''
  Class for a Signal that various representations of a signal (e.g. in time and frequency domains) into one object.
  '''
  # consider making directly iterable/indexable?
  num_signals_made = 0
  def __init__(self, mode='t', data=[0 for _ in range(9)], fs=1, evenlen=True):
  # allow vector of timestamps/'freq'stamps?
    '''
    Initialize Signal given mode, data, and sampling frequency of Signal.

    When mode is 'f', 'w', or 'b', `evenlen` is used to determine the length \
      of the Signal in the time domain (specifically setting `evenlen` to `True`/`False` \
      will result in an even/odd length Signal in the time domain respectively.

    Arguments:
      mode: Domain which the data is originally in. Defaults to 't'. Should be one of \
        't': time (including scaled time) \
        'f': frequency \
        'w': angular frequency \
        'b': normalized angular frequency \
        'beta': same as b \
        'ttilde': scaled time
      data: Raw numerical data
      fs: Sampling rate of data (kHz). Defaults to 1.
      evenlen: Determines whether or not length of Signal in time domain is even if Signal is initialized in frequency domain. Default is True.
    NOTE: fs is not the scaling factor between t/f/w and b/ttilde (that is intrinsic to the filter)
    '''
    self.uid = self.num_signals_made
    Signal.num_signals_made += 1

    data = np.array(data) # helps with MATLAB inputs and brevity of code

    if mode not in ['t', 'f', 'w', 'b', 'beta', 'ttilde']: raise Exception()
    self.mode = mode
    if mode == 'beta':
      self.mode = 'b' # could save a bit of typing
    self.fs = fs

    # having self.length as length of self.mode_t keeps a bit more info (literally)
    if self.mode in ['t', 'ttilde']:
      self.mode_t = data
      self.mode_f = sp.fft.rfft(data)
      self.length = len(data)
    else:
      if self.mode == 'w':
        self.mode_f = data / 2 / np.pi
      else:
        self.mode_f = data
      self.length = 2*len(data)-1-int(evenlen) # -2 if evenlen is True and -1 if evenlen is False
      self.mode_t = sp.fft.irfft(self.mode_f, self.length)

    self.func = (lambda t: self.at_time(t, tolerance=0)) # automatically in time domain, from_function for other domains
    self.timestamps = np.arange(self.length)/fs
    self.freqstamps = scipy.fft.rfftfreq(self.length, 1/fs) # explain why
    if self.mode in ['b', 'beta', 'ttilde']:
      self.freqstamps *= 2 * np.pi

    self.mean = np.mean(self.mode_t)
    self.rms = np.mean([x**2 for x in self.mode_t])**0.5

    self.analytic = sp.signal.hilbert(self.mode_t)
    self.hilbert = np.real(self.analytic)
    self.inst_phase = np.unwrap(np.angle(self.analytic))

  @classmethod
  def from_function(cls, mode='t', func=(lambda x: 0), fs=1, num_samples=9, evenlen=True):
    '''
    Defines Signal from function generating values of Signal. For instance, if `mode` is 't',
    `func` should be a function from time (s) to raw data (usually but not necessarily in kPa).

    When mode is 'f', 'w', or 'b', `evenlen` is used to determine the length \
      of the Signal in the time domain (specifically setting `evenlen` to `True`/`False` \
      will result in an even/odd length Signal in the time domain respectively.

    Arguments:
      mode: Domain data is originally in. See `__init__` for options
      func: Function generating data of Signal
      fs: Sampling rate of data (kHz). Defaults to 1.
      num_samples: Number of places to sample `func` at
      evenlen: Determines whether or not length of Signal in time domain is even if Signal is initialized in frequency domain. Default is True.
    '''
    if mode in ['t', 'ttilde']:
      sample_points = np.arange(num_samples)/fs
    else:
      sample_points = sp.fft.rfftfreq(2*num_samples-1-int(evenlen), 1/fs)
    S = cls(mode=mode, data=func(sample_points), fs=fs, evenlen=evenlen)
    S.func = func
    return S

  @classmethod
  def linear_chirp(cls, f_init=1, w_init=None, f_final=10, w_final=None, fs=1, num_samples=9):
    '''
    Generates linear chirp with initial instantaneous frequency either
    f_init or w_init and final instantaneous frequency either f_final or w_final.

    Only one of {f_init, w_init} and one of {f_final, w_final} needs to be specified.
    An angular frequency overrides the non-angular frequency.

    Arguments:
      f_init: Initial instantaneous frequency
      w_init: Initial angular instantaneous frequency
      f_final: Final instantaneous frequency
      w_final: Final angular instantaneous frequency
      fs: Sampling rate of data. defaults to 1
      num_samples: Number of data points in Signal (in time)
    '''
    fi = f_init
    if w_init:
      fi = w_init / 2 / np.pi

    ff = f_final
    if w_final:
      ff = w_final / 2 / np.pi

    endtime = (num_samples-1)/fs
    return cls.from_function(mode='t', func=(lambda t: np.cos(np.pi*t*(fi*(2-t/endtime) + ff*(t/endtime)))), fs=fs, num_samples=num_samples, evenlen=(num_samples%2==0))
    # evenlen is actually never used in this case since the mode is always 't', but just for completeness

  # @classmethod
  # def from_instantaneous_frequency(cls, freq_func=(lambda x: 0), freqs=None, init_phase=0, fs=1, num_samples=9):
  #   if freqs is None:
  #     # gaussian quadrature?
  #     ws = [2*np.pi*freq_func(i/fs) for i in range(num_samples)]
  #   else:
  #     ws = [2*np.pi*f for f in freqs]

  #   phases = sp.integrate.cumulative_trapezoid(ws, dx=1/fs, initial=0)+init_phase

  #   return cls(mode='t', data=np.cos(phases).tolist(), fs=fs)

  def __iter__(self):
    self.__idx = -1
    return self.get_data(mode=self.mode)

  def __next__(self):
    self.__idx += 1
    try:
      return self.get_data(mode=self.mode)[self.__idx]
    except:
      self.__idx = -1
      raise StopIteration

  def __len__(self):
    return self.length

  def __getitem__(self, idx):
    if idx in ['t', 'f', 'w', 'b', 'beta', 'ttilde']:
      return self.get_data(idx)
    raise Exception('Invalid index into Signal')

  def _pad_time_series_to_same_length(self, other): # should be private
    if self.length >= len(other):
      self_t = self.mode_t
      other_t = np.append(other, [0 for _ in range(self.length-len(other))]).tolist()
    else:
      self_t = np.append(self.mode_t, [0 for _ in range(len(other)-self.length)])
      other_t = other
    return [self_t, other_t]

  def __add__(self, other):
    if isinstance(other, (int, float, np.integer, np.floating)):
      return Signal(mode='t', data=[self.mode_t[i]+other for i in range(self.length)], fs=self.fs)
    else:
      if isinstance(other, Signal):
        if self.fs != other.fs:
          raise Exception('Cannot add two signals with different sampling rates')
        self_t, other_t = self._pad_time_series_to_same_length(other.mode_t)
      else:
        self_t, other_t = self._pad_time_series_to_same_length(other)
      return Signal(mode='t', data=[self_t[i]+other_t[i] for i in range(self.length)], fs=self.fs)

  def __radd__(self, other):
    return Signal(mode='t', data=[other+self.mode_t[i] for i in range(self.length)], fs=self.fs)

  def __sub__(self, other):
    if isinstance(other, (int, float, np.integer, np.floating)):
      return Signal(mode='t', data=[self.mode_t[i]-other for i in range(self.length)], fs=self.fs)
    else:
      if isinstance(other, Signal):
        if self.fs != other.fs:
          raise Exception('Cannot subtract two signals with different sampling rates')
        self_t, other_t = self._pad_time_series_to_same_length(other.mode_t)
      else:
        self_t, other_t = self._pad_time_series_to_same_length(other)
      return Signal(mode='t', data=[self_t[i]-other_t[i] for i in range(self.length)], fs=self.fs)

  def __mul__(self, other):
    if isinstance(other, (int, float, np.integer, np.floating)):
      return Signal(mode='t', data=[self.mode_t[i]*other for i in range(self.length)], fs=self.fs)
    else:
      if isinstance(other, Signal):
        if self.fs != other.fs:
          raise Exception('Cannot multiply two signals with different sampling rates')
        self_t, other_t = self._pad_time_series_to_same_length(other.mode_t)
      else:
        self_t, other_t = self._pad_time_series_to_same_length(other)
      return Signal(mode='t', data=[self_t[i]*other_t[i] for i in range(self.length)], fs=self.fs)

  def __rmul__(self, other):
    return Signal(mode='t', data=[other*self.mode_t[i] for i in range(self.length)], fs=self.fs)

  def __truediv__(self, other):
    if isinstance(other, (int, float, np.integer, np.floating)):
      return Signal(mode='t', data=[self.mode_t[i]/other for i in range(self.length)], fs=self.fs)
    else:
      if isinstance(other, Signal):
        if self.fs != other.fs:
          raise Exception('Cannot divide two signals with different sampling rates')
        self_t, other_t = self._pad_time_series_to_same_length(other.mode_t)
      else:
        self_t, other_t = self._pad_time_series_to_same_length(other)
      return Signal(mode='t', data=[self_t[i]/other_t[i] for i in range(self.length)], fs=self.fs)

  def __mod__(self, other):
    return Signal(mode='t', data=[self.mode_t[i]%other for i in range(self.length)], fs=self.fs)

  def __floordiv__(self, other):
    return Signal(mode='t', data=[self.mode_t[i]//other for i in range(self.length)], fs=self.fs)

  def __pow__(self, other):
    return Signal(mode='t', data=[self.mode_t[i]**other for i in range(self.length)], fs=self.fs)

  def __neg__(self):
    return Signal(mode='t', data=[-self.mode_t[i] for i in range(self.length)], fs=self.fs)

  def __pos__(self):
    return self

  def __abs__(self):
    return Signal(mode='t', data=[abs(self.mode_t[i]) for i in range(self.length)], fs=self.fs)

  def at_time(self, t, tolerance=1e-10):
    '''
    Gets the value of the Signal at time t.

    Fourier interpolates if t is not a multiple of the distance between samples with a tolerance

    Arguments:
      t: time to get value of Signal at (in ms)
      tolerance: if t is within this value of a timestamp, \
        t is clamped to that timestamp and the value at that timestamp is taken
    '''
    num = t * self.fs
    if tolerate(num, eps=tolerance):
      return self.mode_t[round(num)%self.length]

    f = [k/self.length for k in self.mode_f]

    tot = f[0].real
    for idx in range(1, len(f)):
      tot += 2 * (f[idx] * np.exp(2j*np.pi*idx*t/self.length)).real
    if self.length%2 == 0:
      tot -= (f[idx] * np.exp(2j*np.pi*idx*t/self.length)).real
    if tolerate(tot, eps=tolerance):
      tot = round(tot)
    return tot

  def get_data(self, mode='t'):
    '''
    Get the data series in terms of mode (default is in time domain). \
      If Signal was original defined in 't' or 'w' or 'f', the only \
      allowed modes are any of 't'/'w'/'f'. Similarly, Signals originally \
      in 'b'/'beta' or 'ttilde' are only allowed to be accessed in \
      'b'/'beta'/'ttilde' mode

    Arguments:
      mode: Mode to get the data in. Default is 't' (time). See `__init__` for all options.
    '''
    if mode in ['t', 'f', 'w']:
      if self.mode not in ['t', 'f', 'w']:
        raise Exception('Undefined values in requested domain')
      if mode == 't':
        return self.mode_t
      elif mode == 'f':
        return self.mode_f
      else:
        return [2*np.pi*k for k in self.mode_f]
    elif mode in ['b', 'beta', 'ttilde']:
      if self.mode not in ['b', 'beta', 'ttilde']:
        raise Exception('Undefined values in requested domain')
      if mode == 'ttilde':
        return self.mode_t
      else:
        return self.mode_f
    else:
      raise Exception('Invalid domain')

  def resample(self, new_fs, end_time=None):
    '''
    Resamples Signal to new_fs
    '''
    if end_time is None:
      end_time = (self.length-1)*self.fs
    # num_samples = round(end_time/new_fs)+1
    return Signal(data=[self.at_time(i/new_fs) for i in range(round(end_time/new_fs)+1)], fs=new_fs)

  def scale_signal(self, scale_factor, mode=None):
    '''
    Scale signal by scale_factor
    For instance, doubling fs by 2
    '''
    if mode is None:
      if self.mode in ['t', 'f', 'w']:
        mode = 't'
      else:
        mode = 'ttilde'
    elif mode not in ['t', 'ttilde']:
      raise Exception('When scaling a Signal the Signal should end up in t or ttilde mode')
    return Signal(mode=mode, data=self.mode_t, fs=self.fs*scale_factor)

  def envelope_analytic(self):
    '''
    Outputs [upper, lower] where `upper` is the upper envelope \
      of the Signal and `lower` is the lower envelope.
    '''
    centered = self.analytic - self.mean
    inst_amp = abs(centered)
    upper = [x+self.mean for x in inst_amp]
    lower = [self.mean-x for x in inst_amp]
    return [upper, lower]

  def instantaneous_phase(self):
    '''
    Returns instantaneous phase of Signal at `self.timestamps`
    '''
    return self.inst_phase.tolist()

  def instantaneous_freq(self):
    '''
    Returns instantaneous frequency of Signal at `self.timestamps`
    '''
    return np.gradient(self.inst_phase, self.timestamps).tolist()

  def spectral_entropy(self):
    '''
    Returns overall spectral entropy of Signal
    '''
    # really just moving window with window length = self.length
    return self.moving_spectral_entropy(window=self.length)[0]

  def moving_spectral_entropy(self, window_len=9):
    '''
    Returns spectral entropy of a window of Signal as window
    moves over the entire Signal

    Arguments:
      window_len: length of window to calculate spectral entropy
    '''
    Hs = []
    for i in range(self.length-window_len+1):
      data = self.mode_t[i:i+9]
      spectrum = abs(sp.fft.fft(data))**2
      distribution = spectrum/sum(spectrum)
      Hs += [-sum(p*np.log(p) for p in distribution)/np.log(window_len)]
    return Hs

  def spectrogram(self, win=sp.signal.windows.gaussian(30, std=5, sym=True), hop=1, mfft=200, custom_title='Spectrogram', show=True):
    '''
    Generates spectrogram of Signal. Returns [SFFT data, bounds]. \
      Since the window has a small width, the resulting SFFT is \
      actually slightly wider than the original Signal. The proper \
      cut-off bounds to ensure no edge effects from the window are \
      given as the second output.

    Arguments:
      win, hop, mfft: Same as for scipy.signal.ShortTimeFFT.
      custom_title: Optional title of plot. Default is 'Spectrogram'.
      show: `True` if plot is to be shown, `False` otherwise. Default is `True`.
    '''
    N = self.length
    ft = sp.signal.ShortTimeFFT(np.array(win), hop, self.fs, mfft=mfft)
    S = ft.spectrogram(np.array(self.mode_t))
    windowed_S = S[:, ft.lower_border_end[1]:ft.upper_border_begin(N)[1]]
    bounds = (0, ft.delta_t*len(windowed_S[0]), *ft.extent(N)[2:])

    if show:
      fig, ax = plt.subplots()
      img = ax.imshow(abs(windowed_S[::-1]), cmap='viridis', aspect='auto', extent=bounds)
      fig.suptitle('Spectrogram' if custom_title is None else custom_title)
      fig.colorbar(img)
      ax.set_xlabel('Time (ms)')
      if self.mode in ['t', 'w', 'f']:
        ax.set_ylabel('Frequency (kHz)')
      else:
        ax.set_ylabel('Normalized frequency (kHz)')
      plt.show()
    return [windowed_S, bounds]

  def crosscorrelate(self, s2):
    '''
    Crosscorrelates this Signal with other Signal s2
    '''
    if not isinstance(s2, Signal):
      raise Exception('can only cross correlate two signals')
    return np.correlate(self.mode_t, s2.mode_t, mode='full').tolist()

  def __matmul__(self, s2):
    return self.crosscorrelate(s2)

  def autocorrelate(self):
    '''
    Return autocorrelation of Signal (Signal correlated with self).
    '''
    full_corr = np.correlate(self.mode_t, self.mode_t, mode='full')
    half_corr = full_corr[len(full_corr)//2:].tolist()
    return half_corr

  def autocorrelation_plot(self, custom_title='Autocorrelation plot'):
    '''
    Return plot of autocorrelation of Signal

    Arguments:
      custom_title: Optional title of plot. Default is 'Autocorrelation plot'.
    '''
    plt.plot(self.timestamps, self.autocorrelate())
    plt.xlabel('Offset (s)')
    plt.title('Autocorrelation plot' if custom_title is None else custom_title)
    plt.show()

  def plot(self, mode=None, custom_title=None):
    '''
    Plots data series from the indicated mode. If mode is 't'/'ttilde', \
      then the time series is plotted against `self.timestamps`. Otherwise \
      the frequency series is plotted against `self.freqstamps`.

    Similar to `get_data`, Signals initialized in 't'/'w'/'f' can only be plotted \
      in 't'/'w'/'f' and Signals initialized in 'b'/'beta'/'ttilde' can only be \
      plotted in 'b'/'beta'/'ttilde'

    Arguments:
      mode: See `__init__` for modes. Default is 't'/'ttilde' depending on which is defined.
      custom_title: Optional title of plot. NOTE: Default is no title.
    '''
    if mode is None:
      if self.mode in ['t', 'f', 'w']:
        mode = 't'
      else:
        mode = 'ttilde'
    d = self.get_data(mode=mode)
    if mode == 't':
      plt.plot(self.timestamps, d)
      plt.xlabel('Time (s)')
    elif mode == 'ttilde':
      plt.plot(self.timestamps, d)
      plt.xlabel('Normalized time (something)')
    elif mode == 'f':
      plt.plot(scipy.fft.rfftfreq(len(self.mode_t), 1/self.fs), [abs(v) for v in d])
      plt.xlabel('Magnitude of frequency')
    elif mode == 'b' or mode == 'beta':
      plt.plot(scipy.fft.rfftfreq(len(self.mode_t), 1/self.fs), [abs(v) for v in d])
      plt.xlabel('Magnitude of normalized frequency')

    if custom_title is not None:
      plt.title(custom_title)

    plt.show()

  def as_sound(self, filename):
    '''
    Save Signal as sound file

    Attributes:
      filename: Filename to save sound file under.
    '''
    sp.io.wavfile.write(filename, self.fs, np.array(self.mode_t))