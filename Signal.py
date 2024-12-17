import numpy as np
import scipy as sp
import scipy.fft
import matplotlib.pyplot as plt

def tolerate(x, eps=1e-10):
  return abs(x-round(x)) <= eps

class Signal:
  # consider making directly iterable/indexable?
  # make sure spectrogram captions are changed if mode is b/ttilde?
  num_signals_made = 0
  def __init__(self, mode='t', data=[0 for _ in range(9)], fs=1):
  # allow vector of timestamps/'freq'stamps?
    '''
    Signal as object

    Note: lists and Numpy arrays are generally automatically converted to
      Signals with default sample rate 44100

    Arguments:
      mode: domain data is originally in. Defaults to 't'. Should be one of \
        't': time (including scaled time) \
        'f': frequency \
        'w': angular frequency \
        'b': normalized angular frequency \
        'beta': same as b \
        'ttilde': scaled time
      data: raw numerical data
      fs: sampling rate of data. defaults to 1
      **fs is not the scaling factor between t/f/w and b/ttilde (that is intrinsic to the filter)

    Attributes:
      uid: UID of filter

    '''
    self.uid = self.num_signals_made
    Signal.num_signals_made += 1

    data = np.array(data) # helps with MATLAB inputs

    if mode not in ['t', 'f', 'w', 'b', 'beta', 'ttilde']: raise Exception() # get ttilde mode
    self.mode = mode
    if mode == 'beta':
      self.mode = 'b'
    self.fs = fs

    if self.mode in ['t', 'ttilde']:
      self.mode_t = data
      self.mode_f = sp.fft.rfft(data)
    else:
      if self.mode == 'w':
        self.mode_f = data / 2 / np.pi
      else:
        self.mode_f = data
      self.mode_t = sp.fft.irfft(self.mode_f)

    self.func = (lambda t: self.at_time(t, tolerance=0)) # automatically in time domain, from_function for other domains
    self.length = len(self.mode_t) # more info than len(self.mode_f)
    # if timestamps is not None:
    #   self.timestamps = timestamps
    # else:
    #   self.timestamps = [i/fs for i in range(self.length)]
    self.timestamps = np.arange(self.length)/fs
    self.freqstamps = scipy.fft.rfftfreq(self.length, 1/fs) # explain why
    if self.mode in ['b', 'beta', 'ttilde']: # 
      self.freqstamps *= 2 * np.pi

    self.mean = np.mean(self.mode_t)
    self.rms = np.mean([x**2 for x in self.mode_t])**0.5

    self.analytic = sp.signal.hilbert(self.mode_t)
    self.hilbert = np.real(self.analytic)
    self.inst_phase = np.unwrap(np.angle(self.analytic))
    # mag2db?

  @classmethod
  def from_function(cls, mode='t', func=(lambda x: 0), fs=1, num_samples=9):
    '''
    Time is in ms. Frequency is in kHz.

    Arguments:
      mode: domain data is originally in. See __init__ for options
      func: function generating data points of Signal
      fs: sampling rate of data. defaults to 1
      num_samples: number of samples to evaluate func at
    '''
    S = cls(mode=mode, data=[func(x/fs) for x in range(round(num_samples))], fs=fs)
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
      f_init: initial instantaneous frequency
      w_init: initial angular instantaneous frequency
      f_final: final instantaneous frequency
      w_final: final angular instantaneous frequency
      fs: sampling rate of data. defaults to 1
      num_samples: number of data points in Signal (in time)
    '''
    fi = f_init
    if w_init:
      fi = w_init / 2 / np.pi

    ff = f_final
    if w_final:
      ff = w_final / 2 / np.pi

    endtime = (num_samples-1)/fs
    return cls.from_function(mode='t', func=(lambda t: np.cos(np.pi*t*(fi*(2-t/endtime) + ff*(t/endtime)))), fs=fs, num_samples=num_samples) # np.cos(2*np.pi*t*(fi*(1-t/endtime/2) + ff*(t/endtime/2)))

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
    Gets the data at time t.

    Fourier interpolates if t is not a multiple of the distance between samples with a tolerance

    Arguments:
      t: time to get value of Signal at (in ms)
      tolerance: if t is within this value of a timestamp, \
        t is clamped to that timestamp and the value at that timestamp is taken

    # figure these out again
    '''
    num = t * self.fs
    # print(num)
    if tolerate(num, eps=tolerance):
      # print('hi')
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
    Get the data (default is in time domain)

    Arguments:
      mode: what mode to get the data in
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
    Outputs [upper, lower] where *upper* is the upper envelope \
    of the Signal and *lower* is the lower envelope
    '''
    centered = self.analytic - self.mean
    inst_amp = abs(centered)
    upper = [x+self.mean for x in inst_amp]
    lower = [self.mean-x for x in inst_amp]
    return [upper, lower]

  def instantaneous_phase(self):
    '''
    Returns instantaneous phase of Signal at timestamps
    '''
    return self.inst_phase.tolist()

  def instantaneous_freq(self):
    '''
    Returns instantaneous frequency of Signal at timestamps
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
    Generates spectrogram of Signal

    Arguments:
      win: window for spectrogram
      hop: hop
      mfft:
      custom_title: title of spectrogram
      show:
    '''
    N = self.length
    ft = sp.signal.ShortTimeFFT(np.array(win), hop, self.fs, mfft=mfft)
    S = ft.spectrogram(np.array(self.mode_t))
    windowed_S = S[:, ft.lower_border_end[1]:ft.upper_border_begin(N)[1]]
    bounds = (0, ft.delta_t*len(windowed_S[0]), *ft.extent(N)[2:])

    if show:
      fig, ax = plt.subplots()
      img = ax.imshow(abs(windowed_S[::-1]), cmap='viridis', aspect='auto', extent=bounds)
      fig.suptitle(custom_title)
      fig.colorbar(img)
      ax.set_xlabel('Time (s)')
      ax.set_ylabel('Frequency (1/s)')
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

  def autocorrelate(self, custom_title=None, show=True):
    '''
    Return autocorrelation plot
    '''
    full_corr = np.correlate(self.mode_t, self.mode_t, mode='full')
    half_corr = full_corr[len(full_corr)//2:].tolist()
    # if not show:
    #   return half_corr

    if show:
      plt.plot(self.timestamps, half_corr)
      plt.xlabel('Offset (s)')
      plt.title('Autocorrelation plot' if custom_title is None else custom_title)
      plt.show()
    return half_corr

  def plot(self, mode=None, custom_title=None):
    '''
    don't delete this
    '''
    if mode is None:
      if self.mode in ['t', 'f', 'w']:
        mode = 't'
      else:
        mode = 'ttilde'
    d = self.get_data(mode=mode)
    # print('d', d)
    # m = abs(max(d))
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

  # def scatter(self, mode='t'):
  #   '''
  #   don't delete this
  #   '''
  #   d = self.get_data(mode=mode)
  #   plt.scatter([i/self.fs for i in range(len(d))], d, s=3)
  #   plt.show()

  def as_sound(self, filename):
    '''
    Save Signal as sound file
    '''
    sp.io.wavfile.write(filename, self.fs, np.array(self.mode_t))

  # def apply_filter(self, filter, method='default', fs=100):
  #   return filter.solve(self, method=method, fs=fs) # probably not a great sign when importing Filter to typecheck this seems worse

def as_signal(arr) -> Signal:
  if isinstance(arr, Signal):
    return arr
  return Signal(mode='t', data=arr, fs=44100) # should 44.1k be default or 1

if __name__ == "__main__":
  # s1 = Signal(mode='t', data =[1, 4, -2, 3, 0])
  # s2 = Signal(mode='f', data=np.ones(11))
  # s3 = Signal.from_function(mode='t', func=(lambda t: np.sin (5*t**2)), fs=1000, n=2000)
  # s4 = Signal.from_instantaneous_frequency(func=(lambda t: 1+t), init_phase=0, fs=1000, n=5000)
  # s5 = Signal.from_instantaneous_frequency(freqs=[6-i/1000 for i in range(5000)], init_phase=0, fs=1000, n=5000)
  # s1.plot()
  # s2.plot()
  # s3.plot()
  # s4.plot()
  # s5.plot()
  # s6 = Signal.linear_chirp(f_init=1.4, f_final=14, fs=1000, n=2000)
  # s6.plot()

  # s6 = Signal.linear_chirp(f_init=1.5, f_final=9, fs=100, num_samples=200)
  # xaxis = [i/1000 for i in range(len(s6))]
  # plt.plot(xaxis, s6['t'])
  # # plt.plot(xaxis, s6.instantaneous_phase())
  # plt.plot(xaxis[1:-1], s6.instantaneous_freq()[1:-1])
  # plt.show()

  # s6 *= (1+np.sin([i/100 for i in range(2000)]))
  # shift_s6 = s6+1
  # corr = s6 @ shift_s6
  # shift_s6.plot(custom_title='shift_s6')
  # # corr.plot(custom_title='corr')
  # print(len(corr))
  # plt.plot([x/1000 for x in range(-1999, 2000)], corr)
  # plt.title('corr')
  # plt.show()

  # old_s6 = s6
  # new_s6 = old_s6*(1+np.sin([i/100 for i in range(2000)]))
  # new_s6 += 1
  # new_s6.plot()

  # WHY IS THIS ENVELOPE BORKEN

  # xaxis = [i/1000 for i in range(len(new_s6))]
  # plt.plot(xaxis, new_s6['t'])
  # # plt.plot(xaxis, s.inst_amp)
  # upper, lower = new_s6.envelope_analytic()
  # plt.plot(xaxis, upper)
  # plt.plot(xaxis, lower)
  # plt.plot(xaxis, new_s6.instantaneous_phase())
  # plt.plot(xaxis, new_s6.instantaneous_freq())
  # plt.show()

  # s = Signal.from_function(mode='t', func=(lambda t: np.sin(60*t)), fs=200, n=1000)
  # s.plot(custom_title='Original signal')
  # s *= [(1100-i)*(2+np.sin(i/35))/1000 for i in range(1000)]
  # s.plot(custom_title='Modulated signal')
  # xaxis = [i/200 for i in range(1000)]
  # upper, lower = s.envelope_analytic()
  # plt.plot(xaxis, s['t'])
  # plt.plot(xaxis, upper)
  # plt.plot(xaxis, lower)
  # plt.title('Envelope (from analytic signal)')
  # plt.show()

  midC = Signal.linear_chirp(f_init=200, f_final=400, fs=1000, num_samples=1000)
  midC.spectrogram()
  # print(len(midC))
  # print(len(midC.moving_spectral_entropy()))

  # s = Signal.from_function(mode='t', func=(lambda t: np.sin(60*t)), fs=200, n=1000)
  # s.autocorrelate()

  # s = Signal.linear_chirp(f_init=10, f_final=40, fs=2000, n=1000)
  # s.spectrogram(mfft=100)
  # print(s.spectral_entropy())
  # s.plot()
  # # xaxis = [i/200 for i in range(1000)]
  # Hs = s.moving_spectral_entropy()
  # plt.plot(np.linspace(0, 0.5, len(Hs)), Hs)
  # plt.title('Moving spectral entropy')
  # plt.show()

  ###############


  # a = [1, 4, -2, 3, 9, 0, 1, 2, 9, 3, 1]
  # s = Signal(data=a)
  # print(s.mode_f)
  # for i in range(21):
  #   print(i/2, s.at_time(i/2))

  # fs = [(1)*np.exp(1j*i/100) for i in range(100)]
  # s = Signal(mode='f', data=fs)
  # s.plot()

  # s = Signal.from_function(func=(lambda t: np.sin(5*t**2)), fs=100, n=200)
  # # s.plot(mode='f')
  # print(s['f'])
  # print(abs(np.array(s['f'])))
  # s = Signal.linear_chirp(c=2, fs=100, n=200)
  # s.plot()

  # sig = Signal.from_instantaneous_frequency(func=(lambda t: 4-2*t), fs=1000, n=1000)
  # plt.plot(range(len(sig)), sig['t'])

  # sr = 1001
  # # s = Signal.linear_chirp(f_init=5, f_final=0.5, fs=sr, n=10*sr)
  # s = Signal.linear_chirp(f_init=-3.1, f_final=3.1, fs=sr, n=9*sr)
  # s.plot()

  # s *= (1+np.sin([i/sr for i in range(len(s))]))
  # # # s = Signal(mode='t', data=s['t'], fs=1000)

  # xaxis = [i/sr for i in range(len(s))]
  # # plt.plot(xaxis, s['t'])
  # # plt.plot(xaxis, s.inst_amp)
  # # upper, lower = s.envelope_analytic()
  # # plt.plot(xaxis, upper)
  # # plt.plot(xaxis, lower)
  # plt.plot(xaxis, s.instantaneous_phase())
  # plt.show()

  # sr = 1000
  # octave_around_middle_Cish = Signal.linear_chirp(f_init=200, f_final=400, fs=sr, n=1*sr)
  # OC = octave_around_middle_Cish
  # data = OC.spectrogram(hop=5)


  # white_noise = Signal(mode='t', data=np.random.normal(0, 0.5, 22050), fs=44100)
  # white_noise.as_sound('white_noise.wav')

  # sr = 1000
  # num = 1*sr
  # A440 = Signal.linear_chirp(f_init=220, f_final=880, fs=sr, n=num)
  # A440 *= (1.5+np.sin([i/sr for i in range(num)]))
  # A440.plot()
  # # A440 *= Signal.linear_chirp(f_init=220, f_final=880, fs=sr, n=sr)
  # # win = [0.1, 0.5, 1, 0.5, 0.1]
  # # win = [0.1, 1, 1, 1, 0.1]
  # win_num_half = 49
  # win = sp.signal.windows.gaussian(2*win_num_half+1, 1)
  # S = A440.stft(win=win, hop=10)
  # # print(S)
  # lenS = len(S)
  # l = len(S[0])
  # xaxis = [i/sr for i in range(-win_num_half, l-win_num_half)]
  # # moving_average = S[0]
  # # base_freq = S[1]

  # plt.plot([i/sr for i in range(len(A440))], A440['t'])
  # # dot_width = 2
  # # for i in range(0, win_num_half, win_num_half//5):
  # #   fseries = S[i]
  # #   # plt.plot(xaxis, [abs(c) for c in fseries])
  # #   plt.scatter(xaxis, [abs(c) for c in fseries], s=dot_width)
  # plt.show()


  # a = np.correlate(s['t'], s['t'], mode='same')
  # b = np.correlate(s['t'], s['t'], mode='full')
  # c = b[len(b)//2:]
  # plt.plot(range(len(a)), a)
  # plt.plot(range(len(b)), b)
  # c = s.autocorrelate()
  # plt.plot(range(len(c)), c)
  # plt.show()



  # num = 1000
  # func1 = (lambda t: t*(4*(2-t/1) + 2*(t/1)))
  # args1 = [func1(i/num) for i in range(num+1)]
  # plt.plot(range(num+1), args1)
  # func2 = (lambda t: 4-2*t)
  # w0 = func2(0)
  # ws = [2*func2(i/num) for i in range(num+1)]
  # # phases = [0 for _ in range(num+1)]
  # # for i in range(1, num+1):
  # #   phases[i] = phases[i-1] + (ws[i]-ws[i-1])
  # # print(ws)
  # # phases = np.cumsum(ws)/num
  # phases = sp.integrate.cumulative_trapezoid(ws, dx=1/num, initial=0)
  # # phases -= 2*w0/num
  # args2 = phases
  # # args2 = [phases[i]*i/num for i in range(num+1)]
  # plt.plot(range(num+1), args2)
  # plt.show()
  # # print(args1)
  # # print(args2)

