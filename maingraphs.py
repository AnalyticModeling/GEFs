from Filter import *
from FilterBank import *
from Signal import *
from OutputSignals import *
from Cochlea import *
import matplotlib.pyplot as plt
import helpers

### for auditory physicists

def fig3_2019():
  c = Cochlea(Ap=[0.05], bp=[1], Bu=[1.3])
  betas, reals, imags, magns, phases = c.plot_wavenumber(betas=np.geomspace(0.7, 1.2, 10000), show=False)
  fig, (ax1, ax2) = plt.subplots(2, 1)

  ax1.semilogx(betas, reals, ls='--')
  ax1.xaxis.set_major_locator(locator=matplotlib.ticker.LogLocator(subs=(0.7, 1, 1.2)))
  ax1.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
  ax1.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
  ax1.set_ylabel(r'Re{k} (mm⁻¹)')
  ax1.axhline(y=0, color='k', ls=':')
  ax1.axvline(x=1, color='k', ls=':')

  ax2.semilogx(betas, imags, ls='--')
  ax2.xaxis.set_major_locator(locator=matplotlib.ticker.LogLocator(subs=(0.7, 1, 1.2)))
  ax2.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
  ax2.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
  ax2.set_ylabel(r'Im{k} (mm⁻¹)')
  ax2.axhline(y=0, color='k', ls=':')
  ax2.axvline(x=1, color='k', ls=':')
  ax2.set_xlabel('β')
  plt.show()

def fig6_2019():
  c = Cochlea.five_param(types='V', aAp=0.3768, bAp=-0.1366, bp=[0.2, 0.5, 2, 5, 15], aBu=3.714, bBu=0.03123, xs=[i for i in range(5)])
  c.bode_plot(freqs=np.geomspace(0.1, 20, 10000))
  # for fil in c.filters:
  #   print(fil.get_computed_chars())
  # what are the x values?
  # why's the first filter have Bpeak=0

def fig8_2019():
  c = Cochlea.five_param(types='V', aAp=0.3768, bAp=-0.1366, bp=[1, 1, 1, 1, 1], aBu=3.714, bBu=0.03123, xs=[i for i in range(5)])
  pairs = [(1.5, 200), (8, 400), (1.5, 700), (0.3, 400)]
  def tones(t):
    ans = 0
    for i in range(4):
      fi, ti = pairs[i]
      ans += np.exp(-((t-ti)/50)**2) * np.sin(2*np.pi*fi*t)
    return ans
  sig = Signal(mode='t', data=[tones(t/100) for t in range(100000)], fs=100)
  c.signal_response_heatmap(sig)

def fig2_2022():
  fig, axs = plt.subplots(2, 2, constrained_layout=True)

  c1 = Cochlea(Ap=[0.11], bp=[1], Bu=[7], CF0=1)
  betas1, reals1, imags1, magns1, phases1 = c1.plot_wavenumber(betas=np.geomspace(0.5, 2, 10000), show=False)
  c2 = Cochlea(Ap=[0.055], bp=[1], Bu=[7], CF0=10)
  betas2, reals2, imags2, magns2, phases2 = c2.plot_wavenumber(betas=np.geomspace(0.5, 2, 10000), show=False)
  axs[0][0].semilogx(betas1, reals1)
  axs[0][0].semilogx(betas2, reals2, ls='--')
  axs[0][0].xaxis.set_major_locator(locator=matplotlib.ticker.LogLocator(subs=(0.5, 1, 2)))
  axs[0][0].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
  axs[0][0].xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
  axs[0][0].set_ylabel(r'Re{k} (mm⁻¹)')
  axs[0][0].axhline(y=0, color='k', ls=':')
  axs[0][0].axvline(x=1, color='k', ls=':')

  axs[0][1].semilogx(betas1, imags1)
  axs[0][1].semilogx(betas2, imags2, ls='--')
  axs[0][1].xaxis.set_major_locator(locator=matplotlib.ticker.LogLocator(subs=(0.5, 1, 2)))
  axs[0][1].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
  axs[0][1].xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
  axs[0][1].set_ylabel(r'Im{k} (mm⁻¹)')
  axs[0][1].axhline(y=0, color='k', ls=':')
  axs[0][1].axvline(x=1, color='k', ls=':')

  c1 = Cochlea(Ap=[0.11], bp=[1], Bu=[7], CF0=1)
  betas1, reals1, imags1, magns1, phases1 = c1.plot_impedance(betas=np.geomspace(0.5, 2, 10000), show=False)
  c2 = Cochlea(Ap=[0.055], bp=[1], Bu=[7], CF0=10)
  betas2, reals2, imags2, magns2, phases2 = c2.plot_impedance(betas=np.geomspace(0.5, 2, 10000), show=False)
  axs[1][0].semilogx(betas1, reals1)
  axs[1][0].semilogx(betas2, reals2, ls='--')
  axs[1][0].xaxis.set_major_locator(locator=matplotlib.ticker.LogLocator(subs=(0.5, 1, 2)))
  axs[1][0].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
  axs[1][0].xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
  axs[1][0].set_ylabel(r'Re{Z}/(2πρl)')
  axs[1][0].axhline(y=0, color='k', ls=':')
  axs[1][0].axvline(x=1, color='k', ls=':')
  axs[1][1].set_xlabel('β')

  axs[1][1].semilogx(betas1, imags1)
  axs[1][1].semilogx(betas2, imags2, ls='--')
  axs[1][1].xaxis.set_major_locator(locator=matplotlib.ticker.LogLocator(subs=(0.5, 1, 2)))
  axs[1][1].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
  axs[1][1].xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
  axs[1][1].set_ylabel(r'Im{Z}/(2πρl)')
  axs[1][1].axhline(y=0, color='k', ls=':')
  axs[1][1].axvline(x=1, color='k', ls=':')
  axs[1][1].set_xlabel('β')
  plt.show()

def fig1_2024():
  Ap = 0.055
  bp = 1
  Bu = 7
  fil = Filter(Ap=Ap, bp=bp, Bu=Bu)
  P = fil.filter.tf
  Psharp = fil.filter.sharp_approximation()[0]
  p = 1j*bp - Ap
  k = (lambda beta: Bu * (1/(1j*beta - p) + 1/(1j*beta - p.conjugate())))
  ksharp = (lambda beta: Bu * (1/(1j*beta - p)))

  freqs = np.geomspace(0.1, 2, 10000)
  responses = [P(1j*freqs), Psharp(1j*freqs), k(freqs), ksharp(freqs)]

  fig, axs = plt.subplots(2, 2, constrained_layout=True)
  Pdb = helpers.mag2db(abs(responses[0]))
  Psharpdb = helpers.mag2db(abs(responses[1]))
  axs[0][0].semilogx(freqs, Pdb-max(Pdb))
  axs[0][0].semilogx(freqs, Psharpdb-max(Psharpdb), ls='--')
  axs[0][0].set_title(r'mag{P} (dB)')

  Pcyc = np.unwrap(np.angle(responses[0])) / (2 * np.pi)
  Psharpcyc = np.unwrap(np.angle(responses[1])) / (2 * np.pi)
  axs[1][0].semilogx(freqs, Pcyc-Pcyc[0])
  axs[1][0].semilogx(freqs, Psharpcyc-Psharpcyc[0], ls='--')
  axs[1][0].set_title(r'phase{P} (cyc)')
  axs[1][0].set_xlabel('β')

  kre = [z.real for z in responses[2]]
  ksharpre = [z.real for z in responses[3]]
  axs[0][1].semilogx(freqs, kre)
  axs[0][1].semilogx(freqs, ksharpre, ls='--')
  axs[0][1].set_title(r'Re{kᵦ}')

  kim = [z.imag for z in responses[2]]
  ksharpim = [z.imag for z in responses[3]]
  axs[1][1].semilogx(freqs, kim)
  axs[1][1].semilogx(freqs, ksharpim, ls='--')
  axs[1][1].set_title(r'Im{kᵦ}')
  axs[1][1].set_xlabel('β')
  plt.show()

def fig3_2024():
  filP = Filter(Bpeak=1, Nbeta=19.1, Qerb=25.9)
  filV = Filter(type='V', Bpeak=1, Nbeta=19.1, Qerb=25.9)
  print(filP.get_computed_chars())
  filP.characteristic_error()

  P = filP.filter.tf
  Psharp = filP.filter.sharp_approximation()[0]
  V = filV.filter.tf

  freqs = np.geomspace(0.1, 8, 10000)
  responses = [P(1j*freqs), Psharp(1j*freqs), V(1j*freqs)]

  fig, axs = plt.subplots(1, 2, constrained_layout=True)
  Pdb = helpers.mag2db(abs(responses[0]))
  Psharpdb = helpers.mag2db(abs(responses[1]))
  Vdb = helpers.mag2db(abs(responses[2]))
  axs[0].semilogx(freqs, Pdb-max(Pdb))
  axs[0].semilogx(freqs, Psharpdb-max(Psharpdb), ls='--')
  axs[0].semilogx(freqs, Vdb-max(Vdb), ls=':')
  axs[0].set_ylabel(r'magnitude (dB)')
  axs[0].set_xlabel('β')

  Pcyc = np.unwrap(np.angle(responses[0])) / (2 * np.pi)
  Psharpcyc = np.unwrap(np.angle(responses[1])) / (2 * np.pi)
  Vcyc = np.unwrap(np.angle(responses[2])) / (2 * np.pi)
  axs[1].semilogx(freqs, Pcyc-Pcyc[0])
  axs[1].semilogx(freqs, Psharpcyc-Psharpcyc[0], ls='--')
  axs[1].semilogx(freqs, Vcyc-Vcyc[0], ls=':')
  axs[1].set_ylabel(r'phase (cyc)')
  axs[1].set_xlabel('β')
  plt.show()

def fig5_2024():
  f1 = Filter(Bpeak=1, Qn=7.8, Sbeta=1.7e3, n=3)
  f2 = Filter(Bpeak=1.5, Qn=12, Sbeta=1.7e3, n=3)
  f3 = Filter(Bpeak=2, Qn=16, Sbeta=1.7e3, n=3)
  fall = Filter.multiband_chars(Bpeak=[1, 1.5, 2], Qn=[7.8, 12, 16], Sbeta=1.7e3, n=3)
  betas = np.linspace(0.1, 2.5, 10000)

  minidx1 = 4791
  minidx2 = 6874
  data = [fil.bode_plot(freqs=betas, show=False) for fil in [f1, f2, f3, fall]]

  fig, axs = plt.subplots(2, 2, constrained_layout=True)

  for d in data:
    axs[0][0].plot(d[0], d[1])
    axs[0][1].plot(d[0], d[2])
  gs = axs[1][0].get_gridspec()
  for ax in axs[1]:
    ax.remove()
  axbig = fig.add_subplot(gs[1, :])

  c1 = helpers.computedfiltercharacteristics(tfunc=fall.filter.tf, betas=betas[:minidx1], n=3)
  c2 = helpers.computedfiltercharacteristics(tfunc=fall.filter.tf, betas=betas[minidx1:minidx2], n=3)
  c3 = helpers.computedfiltercharacteristics(tfunc=fall.filter.tf, betas=betas[minidx2:], n=3)

  labels = [r'1:$\beta_{peak}$', r'2:$\beta_{peak}$', r'3:$\beta_{peak}$', '1:Q$_3$', '2:Q$_3$', '3:Q$_3$', r'1:S$_\beta$', r'2:S$_\beta$', r'3:S$_\beta$']
  expected = [1, 1.5, 2, 7.8, 12, 16, 1.7e3, 1.7e3, 1.7e3]
  estimated = [c1['Bpeak'], c2['Bpeak'], c3['Bpeak'], c1['Qn'], c2['Qn'], c3['Qn'], c1['Sbeta'], c2['Sbeta'], c3['Sbeta']]
  errors = [abs(estimated[i]/expected[i]-1) for i in range(9)]

  axbig.bar(labels, errors)
  plt.show()

def fig1_pending():
  fils = [Filter(Ap=0.05, bp=1, Bu=v) for v in [2, 2.5, 3]]
  freqs = np.geomspace(0.8, 1.2, 10000)
  responses = [fil.filter.tf(1j*freqs) for fil in fils]

  fig, axs = plt.subplots(1, 2, constrained_layout=True)
  twodb = helpers.mag2db(abs(responses[0]))
  twopointfivedb = helpers.mag2db(abs(responses[1]))
  threedb = helpers.mag2db(abs(responses[2]))
  axs[0].semilogx(freqs, twodb-max(twodb))
  axs[0].semilogx(freqs, twopointfivedb-max(twopointfivedb), ls='--')
  axs[0].semilogx(freqs, threedb-max(threedb), ls=':')
  axs[0].set_ylabel('Magnitude (dB)')
  axs[0].set_xlabel(r'$\beta$')

  twocyc = np.unwrap(np.angle(responses[0])) / (2 * np.pi)
  twopointfivecyc = np.unwrap(np.angle(responses[1])) / (2 * np.pi)
  threecyc = np.unwrap(np.angle(responses[2])) / (2 * np.pi)
  axs[1].semilogx(freqs, twocyc-twocyc[0])
  axs[1].semilogx(freqs, twopointfivecyc-twopointfivecyc[0], ls='--')
  axs[1].semilogx(freqs, threecyc-threecyc[0], ls=':')
  axs[1].set_ylabel('Phase (cyc)')
  axs[1].set_xlabel(r'$\beta$')
  plt.show()

def fig2_pending():
  Ncycs = []
  Qerbs = []
  Q10s = []
  Q3s = []
  Q15s = []
  for Bu in range(2, 9):
    fil = Filter(Ap=0.05, bp=1, Bu=Bu)
    Ncycs += [fil.get_computed_chars()['Nbeta']]
    Qerbs += [fil.get_computed_chars()['Qerb']]
    Q10s += [fil.get_computed_chars()['Qn']]
    Q3s += [fil.get_computed_chars()['Qn2']]
    Q15s += [helpers.computedfiltercharacteristics(fil.filter.tf, n=15)['Qn']]
  plt.plot(range(2, 9), Ncycs, label=r'N$_{cyc}$')
  plt.plot(range(2, 9), Qerbs, label=r'Q$_{erb}$')
  plt.plot(range(2, 9), Q10s, label=r'Q$_{10}$')
  plt.plot(range(2, 9), Q3s, label=r'Q$_{3}$')
  plt.plot(range(2, 9), Q15s, label=r'Q$_{15}$')
  plt.axhline(y=0, color='k', ls=':')
  plt.legend()
  plt.show()

def fig3_pending():
  fil = Filter(Ap=0.05, bp=1, Bu=2, cf=1)
  pairs = [(1, 20), (5, 50), (7/8, 70), (1/5, 40)]
  def tones(t):
    ans = 0
    for i in range(4):
      fi, ti = pairs[i]
      ans += np.exp(-((t-ti)/5)**2) * np.sin(2*np.pi*fi*t)
    return ans
  fs = 100
  sig = Signal(mode='t', data=[tones(t/fs) for t in range(fs*100)], fs=fs)
  outsig = fil.solve(sig, method='tf')

  fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True)
  ax1.plot(sig.timestamps, sig.mode_t)
  ax1.set_ylabel('input')
  ax2.plot(outsig.timestamps, outsig.mode_t)
  ax2.set_ylabel('output')
  fig.suptitle('Processing using non-integer TF representation, B$_u$=2.5e+00, A$_p$=5.0e-02, b$_p$=1')
  plt.show()

def fig4_pending():
  fils1 = [Filter(Ap=Ap, bp=1, Bu=Bu) for Ap, Bu in [(0.15, 3), (0.15, 5), (0.045, 5)]]
  fils2 = [Filter(Ap=0.15, bp=1, Bu=Bu) for Bu in [1.5, 2, 2.5, 3]]

  fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True)

  for fil in fils1:
    timestamps = [i/10 for i in range(2000)]
    _, ir = fil.impulse_response_plot(times=timestamps, show=False)
    maxir = max(ir)
    ax1.plot(timestamps, [v/maxir for v in ir])
  ax1.set_title('dependence of behavior of h on values of constants for integer Bᵤ')
  ax1.set_xlabel(r'h($\widetilde{t}$)') # \widetilde?
  ax1.set_ylabel(r'$\widetilde{t}$')
  ax1.axhline(y=0, color='k', ls=':')

  for fil in fils2:
    timestamps = [i/10 for i in range(1000)]
    _, ir = fil.impulse_response_plot(times=timestamps, show=False)
    maxir = max(ir)
    ax2.plot(timestamps, [v/maxir for v in ir])
  ax2.set_title('h and phase of oscillatory component for integer and half-integer Bᵤ')
  ax2.set_xlabel('h($\widetilde{t}$)')
  ax2.set_ylabel('$\widetilde{t}$')
  ax2.axhline(y=0, color='k', ls=':')
  plt.show()

def fig7_pending():
  fil = Filter(Ap=0.1, bp=1, Bu=1.75, cf=1)
  sig = Signal.linear_chirp(f_init=-2, f_final=2, fs=30, num_samples=3000) # set this to be able to be ttilde
  sig = Signal(mode='ttilde', data=sig.get_data('t'), fs=30)
  # sig.plot()
  # sig.plot(mode='beta')
  sol = fil.solve(sig, method='tf')
  # sol.plot()
  # sig.spectrogram()
  # sol.spectrogram()

  fig, axs = plt.subplots(2, 2, constrained_layout=True)
  axs[0][0].plot(sig.timestamps, sig['ttilde'])
  axs[0][1].plot(sol.timestamps, sol['ttilde'])
  winsig, sigbound = sig.spectrogram(show=False)
  winsol, solbound = sol.spectrogram(show=False)
  axs[1][0].imshow(abs(winsig[0:21][::-1]), cmap='gray', aspect='auto', extent=(0, sigbound[1], 0, sigbound[3]/5))
  axs[1][1].imshow(abs(winsol[0:21][::-1]), cmap='gray', aspect='auto', extent=(0, solbound[1], 0, solbound[3]/5))
  plt.show()

def fig8_pending():
  sig = Signal.from_function(mode='ttilde', func=(lambda t: t*np.cos(10*t)*np.exp(-t/2) + t**3*np.exp(-t)*np.cos(t)), fs=10, num_samples=1000)
  pairs = [(1, 20), (5, 50), (7/8, 70), (1/5, 40)]
  def tones(t):
    ans = 0
    for i in range(4):
      fi, ti = pairs[i]
      ans += np.exp(-((t-ti)/5)**2) * np.sin(2*np.pi*fi*t)
    return ans
  fs = 10
  sig = Signal(mode='ttilde', data=[tones(t/fs) for t in range(fs*100)], fs=fs)

  fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True)
  ax1.plot(sig.timestamps, sig.get_data('ttilde'))
  ax1.set_ylabel('input')

  fil = Filter(Ap=0.15, bp=1, Bu=5)
  for method in ['tf', 'ir', 'ode', 'fde']:
    sol = fil.solve(sig, method=method).get_data('ttilde')
    maxsol = max(sol)
    ax2.plot(sig.timestamps, [v/maxsol for v in sol], ls=':', label=method) # add legend
  ax2.set_ylabel('output')
  ax2.set_ylabel('t˜')
  fig.suptitle('equivalence of representations Aₚ=0.15, bₚ=1, Bᵤ=5')
  plt.show()

def fig9_pending():
  Ap = 0.04
  bp = 1
  Bu = 1.5
  sig = Signal.from_function(mode='ttilde', func=(lambda t: np.exp(-Ap*t)*scipy.special.jn(0, t*bp)), fs=10, num_samples=1000)
  fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True)
  ax1.plot(sig.timestamps, sig.get_data('ttilde'))
  ax1.set_ylabel('input')

  fil = Filter(Ap=Ap, bp=bp, Bu=Bu)
  for method in ['tf', 'ir', 'fde']:
    sol = fil.solve(sig, method=method).get_data('ttilde')
    maxsol = max(sol)
    ax2.plot(sig.timestamps, [v/maxsol for v in sol], ls=':')
  ax2.set_ylabel('output')
  ax2.set_ylabel('t˜')
  fig.suptitle('equivalence of representations Aₚ=0.04, bₚ=1, Bᵤ=1.5')
  plt.show()

if __name__ == "__main__":
  fig3_2019()
  fig6_2019()
  fig8_2019()
  fig2_2022()
  fig1_2024()
  fig3_2024()
  fig5_2024()
  fig1_pending()
  fig2_pending()
  fig3_pending()
  fig4_pending()
  fig7_pending()
  fig8_pending()
  fig9_pending()

  # example for all four plots for cochlea
  # one example of filter using rational characteristics
  # design single filter using Qerb and N and whatever and calculate desired errors
