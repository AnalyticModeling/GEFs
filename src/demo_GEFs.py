import numpy as np
from matplotlib import pyplot as plt

from GEFs_core.Filter import Filter
from GEFs_core.FilterBank import FilterBank
from GEFs_core.Signal import Signal


# import scipy.signal as spsig

def filter_init():
  f1 = Filter(tf=(lambda s: 1/(1+s+s**2)))
  f1.bode_plot()
  f2 = Filter(coeffs=[[1], [1, 1, 1]])
  f2.bode_plot()
  f3 = Filter(roots=[[1], [1+2j, 1-2j]])
  f3.bode_plot()
  f4 = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5)
  f4.bode_plot()
  f5 = Filter(Bpeak=1.5, Nf=1.11, phiaccum=3.5, cf=10)
  f5.bode_plot()
  f6 = Filter(Ap=0.1, bp=1, Bu=3)
  f6.bode_plot()

def filter_multiband_consts():
  f = Filter.multiband_consts(Ap=0.1, bp=[0.5, 1, 1.5], Bu=[3, 5, 7], peak_magndb=1)
  f.bode_plot()

def filter_multiband_chars():
  f = Filter.multiband_chars(Bpeak=[0.5, 1, 1.5], Nbeta=[15, 10, 5], phiaccum=3.5)
  f.bode_plot()

def filter_get_computed_chars():
  f = Filter(Bpeak=1, Nbeta=11.1, phiaccum=3.5)
  print(f.get_computed_chars())

def filter_get_computed_unnormalized_chars():
  f = Filter(Bpeak=1, Nf=11.1, phiaccum=3.5, cf=1)
  print(f.get_computed_unnormalized_chars())

def filter_get_orig_chars():
  f = Filter(Bpeak=1, Nbeta=11.1, phiaccum=3.5)
  print(f.get_orig_chars())

def filter_get_consts():
  # get from paper
  f = Filter(Ap=0.1, bp=1, Bu=3)
  print(f.get_consts())

def filter_solve():
  # f = Filter(type='P', Bpeak=1, Nbeta=11.1, phiaccum=3.5)
  # Define the filter
  f = Filter(Ap=0.1, bp=1.0, Bu=2, cf=1)

  # Define the input signal function
  func = lambda t: np.exp(-1 / 15 * (t - 20) ** 2) * np.cos(t)

  # Create the signal
  sig = Signal(mode='ttilde', data=[func(t/20) for t in range(2000)], fs=20)

  # Solve using different methods
  anstf = f.solve(sig, method='tf')
  anstf /= max(anstf.mode_t)
  # anstf.plot(custom_title='tf solve')

  ansir = f.solve(sig, method='ir')
  ansir /= max(ansir.mode_t)
  # ansir.plot(custom_title='ir solve')

  ansode = f.solve(sig, method='ode')
  ansode /= max(ansode.mode_t)
  # ansode.plot(custom_title='ode solve')

  ansfde = f.solve(sig, method='fde')
  ansfde /= max(ansfde.mode_t)
  # ans.plot(custom_title='fde solve')

  # Create a figure with two subplots side-by-side
  plt.figure(figsize=(12, 5))
  plt.suptitle("Comparison of filter output using different methods", fontsize=14)

  # Plot 1: Input Signal
  plt.subplot(1, 2, 1)
  plt.plot(sig.timestamps, sig.mode_t)
  plt.xlabel("Normalized time (ms/kHz)")
  plt.ylabel("Input Signal")

  # Plot 2: Output Methods Comparison
  plt.subplot(1, 2, 2)
  plt.plot(anstf.timestamps, anstf.mode_t, label='tf')
  plt.plot(ansir.timestamps, ansir.mode_t, ls='--', label='ir')
  plt.plot(ansode.timestamps, ansode.mode_t, label='ode')
  plt.plot(ansfde.timestamps, ansfde.mode_t, ls='--', label='fde')

  plt.xlabel("Normalized time (ms/kHz)")
  plt.ylabel("Output Signal")
  plt.legend()
  plt.show()

# def filter_solve_t_vs_ttilde():
#   f = Filter(Ap=0.01, bp=1.0, Bu=3, cf=1)
#   func = lambda t: np.exp(-1/15*(t-20)**2) * np.cos(t)
#   timestamps = [i/10 for i in range(1000)]
#   sig1 = Signal(mode='t', data=[func(t) for t in timestamps], fs=10)
#   sig2 = Signal(mode='ttilde', data=[func(t) for t in timestamps], fs=10/2/np.pi)
#   ans1 = f.solve(sig1, method='tf')
#   ans2 = f.solve(sig2, method='tf')

#   # the following two signals will look different when plotted togther because one is plotted in t and the other in ttilde
#   plt.plot(ans1.timestamps, ans1.mode_t, ls='--', label='t')
#   plt.plot(ans2.timestamps, ans2.mode_t, ls=':', label='ttilde')
#   plt.legend()
#   plt.show()

def filter_bode():
  f1 = Filter(tf=(lambda s: 1/(1+s+s**2)))
  f1.bode_plot()
  f2 = Filter(coeffs=[[1], [1, 1/2, 1/4]])
  f2.bode_plot()
  f3 = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5)
  f3.bode_plot()
  # f4 = Filter(fpeak=1.5, Nbeta=11.1, phiaccum=3.5, cf=1.5)
  # f4.bode_plot(freqs=)

def filter_frequency_real_imag():
  f1 = Filter(tf=(lambda s: 1/(1+s+s**2)))
  f1.frequency_real_imag_plot()
  f2 = Filter(coeffs=[[1], [1, 1/2, 1/4]])
  f2.frequency_real_imag_plot()
  f3 = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5)
  f3.frequency_real_imag_plot()

def filter_nichols():
  f1 = Filter(tf=(lambda s: 1/(1+s+s**2)))
  f1.nichols_plot()
  f2 = Filter(coeffs=[[1], [1, 1/2, 1/4]])
  f2.nichols_plot()
  f3 = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5)
  f3.nichols_plot()

def filter_nyquist():
  f1 = Filter(tf=(lambda s: 1/(1+s+s**2)))
  f1.nyquist_plot()
  f2 = Filter(coeffs=[[1], [1, 1/2, 1/4]])
  f2.nyquist_plot()
  f3 = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5)
  f3.nyquist_plot()

def filter_ir():
  f1 = Filter(tf=(lambda s: 1/(1+s+s**2)))
  f1.impulse_response_plot(times=np.linspace(0, 10, 300),custom_title="$H(s) = 1/(1+s+s^2)$")

  f2 = Filter(coeffs=[[1], [1, 1/2, 1/4]])
  f2.impulse_response_plot(times=np.linspace(0, 10, 300),custom_title="$H(s) = 1 / [1 + (1/2)s + (1/4)s^2]$")

  f3 = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5)
  f3.impulse_response_plot(times=np.linspace(0, 100, 1000),custom_title="Bpeak=1.5, Nbeta=11.1, phiaccum=3.5")

  f4 = Filter(Ap=0.01, bp=1.0, Bu=3, cf=1/2/np.pi)
  f4.impulse_response_plot(times=np.linspace(0, 100, 1000),custom_title="Ap=0.01, bp=1.0, Bu=3, cf=1/2/np.pi")
  # need to make unslow

def filter_pz():
  f1 = Filter(type='V', Bpeak=1.5, Nbeta=11.1, phiaccum=3.5)
  f1.pole_zero_plot()
  f2 = Filter(coeffs=[[1, 2], [1, 1/2, 1/4]])
  f2.pole_zero_plot()
  f3 = Filter(tf=(lambda s: 1/(1+s+s**2)))
  # f3.pole_zero_plot() -> Error

def filter_Qns():
  fil = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5)
  fil.Qn_plot()

def filter_characteristic_error():
  fil = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5)
  fil.characteristic_error()
  # print(fil.characteristic_error())

def signal_init():
  s1 = Signal(mode='t', data=[(1+i/100)*np.sin(i/10) for i in range(100)])
  s1.plot(custom_title="y(t) = (1 + t/100) * sin(t/10)")
  s2 = Signal(mode='f', data=[i/10 for i in range(20)])
  s2.plot(custom_title="Y(f) = f / 10")
  s3 = Signal.from_function(mode='t', func=(lambda x: (1+x/100)*np.sin(x/10)), num_samples=100)
  s3.plot(custom_title="y(t) = (1 + t/100) * sin(t/10)")
  s4 = Signal.from_function(mode='f', func=(lambda x: x/10), num_samples=20)
  s4.plot(custom_title="Y(f) = f / 10")
  s5 = Signal.linear_chirp(f_init=1.5, f_final=6, fs=100, num_samples=800)
  s5.plot(custom_title=" Linear Chirp [1.5Hz to 6Hz]")
  s6 = Signal.from_instantaneous_frequency(freq_func=(lambda x: x/10), fs=20, num_samples=200)
  s6.plot(custom_title="Instantaneous Frequency: f(t) = t/10")

def signal_arith():
  sig = Signal.linear_chirp(f_init=1.5, f_final=6, fs=100, num_samples=800)
  sig.plot(custom_title="Original Signal: y(t) Chirp (1.5Hz to 6Hz)")
  (sig+1).plot(custom_title="y(t) + 1")
  (sig*sig).plot(custom_title="y(t)²")

def signal_at_time():
  sig = Signal.from_function(func=(lambda x: np.sin(x/1.3)), fs=10, num_samples=100)
  sig.plot(custom_title="y(t) = sin(t / 1.3)")
  print(f"Value at 8.1s: {sig.at_time(8.1)}")
  print(f"Value at 2.6πs: {sig.at_time(2.6*np.pi)}")
  print(f"Value at 8.2s: {sig.at_time(8.2)}")

def signal_get_data():
  sig = Signal.linear_chirp(f_init=1.5, f_final=6, fs=100, num_samples=800)
  print(sig.get_data('t'))
  print(sig.get_data('f'))

def signal_resample():
  sig = Signal.linear_chirp(f_init=1.5, f_final=6, fs=100, num_samples=800)
  sig.plot(custom_title="Original Signal: y(t) Chirp (fs = 100 Hz) ")
  sig.resample(new_fs=77).plot(custom_title="Resampled Signal: Chirp (fs = 77 Hz)")

def signal_envelope_analytic():
  sig = Signal.linear_chirp(f_init=1.5, f_final=6, fs=100, num_samples=800)
  new_sig = 1 + sig*(1.5+np.sin([i/65 for i in range(800)]))
  xaxis = [i/100 for i in range(len(new_sig))]
  upper, lower = new_sig.envelope_analytic()
  plt.plot(xaxis, new_sig['t'])
  plt.plot(xaxis, upper)
  plt.plot(xaxis, lower)
  plt.xlabel('Time (ms)')
  plt.ylabel('Linear Chirp Signal (1.5Hz to 6Hz)')
  plt.title("Analytic Envelope Detection: y(t) = 1 + Chirp * (1.5 + sin(t/65))")
  plt.show()

def signal_instantaneous_phase():
  s = Signal.linear_chirp(f_init=200, f_final=400, fs=200, num_samples=1000)
  xaxis = [i/200 for i in range(1000)]
  phases = s.instantaneous_phase()
  plt.plot(xaxis, phases)
  plt.xlabel('Time (s)')
  plt.ylabel('Phase (rad)')
  plt.title('Instantaneous Phase of Linear Chirp')
  plt.show()


def filterbank_add():
  fs = FilterBank(topology='parallel', Ap=[0.1, 0.1], bp=1, Bu=[3, 3], cf=[0.5, 1.0])
  fs.bode_plot()
  fs.add(Filter(Ap=0.1, bp=1.5, Bu=3), source=fs.filters[-1])
  fs.bode_plot()

def filterbank_process_signal():
  # 1. Setup FilterBank
  fb = FilterBank(topology='parallel', Ap=[0.1, 0.1, 0.1], bp=1, Bu=[3, 3, 3], cf=[0.5, 1.0, 1.5])

  # 2. Process the Signal
  chirp = Signal.linear_chirp(f_init=1, f_final=10, fs=100, num_samples=200)
  os = fb.process_signal(chirp)

  # Output for the first filter
  os.outsignals[0].plot(custom_title="FilterBank Output 1 (cf=0.5 kHz)")

  # Output for the second filter
  os.outsignals[1].plot(custom_title="FilterBank Output 2 (cf=1 kHz)")

  # Output for the third filter
  os.outsignals[2].plot(custom_title="FilterBank Output 3 (cf=1.5 kHz)")

  # or this way :
  # for sig in os.outsignals:
  # sig.plot()

def filterbank_bode():
  fs = FilterBank(topology='parallel', Ap=[0.1, 0.05, 0.01], bp=1, Bu=[3, 3, 3], cf=[0.5, 1.0, 1.5])
  fs.bode_plot()

  # fs = FilterBank(topology='parallel', Ap=[0.1, 0.1, 0.1], Bu=[3, 3, 3], fpeak=[0.5, 1.0, 1.5])

  FilterBank(topology='parallel', Ap=[0.1, 0.05, 0.01], bp=[0.5, 1.0, 1.5], Bu=[3, 3, 3]).bode_plot()

def outputsignals_init_kindof():
  fs = FilterBank(topology='parallel', Ap=[0.1, 0.1, 0.1], bp=[0.5, 1.0, 1.5], Bu=[3, 3, 3])
  # Input Signal
  insig = Signal.linear_chirp(f_init=1, f_final=10, fs=100, num_samples=200)
  insig.plot(custom_title="Input Signal: Linear Chirp (1Hz to 10Hz)")
  os = fs.process_signal(insig)

  os.outsignals[0].plot(custom_title="FilterBank Output 1 (Ap=0.1, bp=0.5, Bu=3)")

  os.outsignals[1].plot(custom_title="FilterBank Output 2 (Ap=0.1, bp=1.0, Bu=3)")

  os.outsignals[2].plot(custom_title="FilterBank Output 3 (Ap=0.1, bp=1.5, Bu=3)")

  # or this way :
  # for sig in os.outsignals:
  # sig.plot()

# def outputsignals_readfile():
#   example = 'red_truth'
#   samplerate, data = sp.io.wavfile.read(f'test_signals/{example}.wav')
#   resample_factor = data.shape[0]//1000
#   samplerate /= resample_factor
#   samplerate /= 1000
#   if len(data.shape) == 1:
#     mono = data
#   else:
#     mono = [sum(d) for d in data]
#   mono = [mono[i] for i in range(0, len(mono), resample_factor)]
#   print('len:', len(mono))
#   # mono = mono[:10000]

#   fs = FilterBank(topology='parallel', Ap=[0.1, 0.1, 0.1], bp=[0.5, 1.0, 1.5], Bu=[3, 3, 3])
#   insig = Signal(mode='t', data=mono, fs=samplerate)
#   insig.plot()
#   insig.spectrogram()
#   os = fs.process_signal(insig)
#   for sig in os.outsignals:
#     sig.plot()
#   os.correlogram()

def outputsignal_autocorrelates():
  fs = FilterBank(topology='parallel', Ap=[0.1, 0.1, 0.1], bp=[0.5, 1.0, 1.5], Bu=[3, 3, 3])
  os = fs.process_signal(Signal.linear_chirp(f_init=1, f_final=10, fs=100, num_samples=200))
  os.autocorrelates()

def outputsignal_correlate_with():
  fs = FilterBank(topology='parallel', Ap=[0.1, 0.1, 0.1], bp=[0.5, 1.0, 1.5], Bu=[3, 3, 3])
  os = fs.process_signal(Signal.linear_chirp(f_init=1, f_final=10, fs=100, num_samples=200))
  sig = Signal.linear_chirp(f_init=2, f_final=5, fs=100, num_samples=200)
  os.correlate_with(sig)

def outputsignals_correlogram():
  fs = FilterBank(topology='parallel', Ap=[0.1, 0.1, 0.1], bp=[0.5, 1.0, 1.5], Bu=[3, 3, 3])
  os = fs.process_signal(Signal.linear_chirp(f_init=1, f_final=10, fs=100, num_samples=200))
  os.correlogram()

if __name__ == "__main__":

   filter_init()
  # filter_multiband_consts()
  # filter_multiband_chars()
  # filter_get_computed_chars()
  # filter_get_computed_unnormalized_chars()
  # filter_get_orig_chars()
  # filter_get_consts()
  # filter_solve()
  # filter_bode()
  # filter_frequency_real_imag()
  # filter_nichols()
  # filter_nyquist()
  # filter_ir()
  # filter_pz()
  # filter_Qns()
  # filter_characteristic_error()
  # signal_init()
  # signal_arith()
  # signal_at_time()
  # signal_get_data()
  # signal_resample()
  # signal_envelope_analytic()
  # signal_instantaneous_phase()
  # filterbank_add()
  # filterbank_process_signal()
  # filterbank_bode()
  # outputsignals_init_kindof()
  # outputsignal_autocorrelates()
  # outputsignal_correlate_with()
  # outputsignals_correlogram()