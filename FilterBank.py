from Filter import Filter
from Signal import Signal
from OutputSignals import OutputSignals
from RootedTree import RootedTree

import warnings
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker

class FilterBank:
  def __init__(self, topology=None, filters=None, type='P', ir=None, tf=None, coeffs=None, roots=None, Ap=None, bp=None, Bu=None, gain_const=None, peak_magndb=None, Bpeak=None, fpeak=None, phiaccum=None, Nbeta=None, Nf=None, Qerb=None, ERBbeta=None, ERBf=None, Qn=None, Qn2=None, BWndBbeta=None, BWndBf=None, BWn2dBbeta=None, BWn2dBf=None, Sbeta=None, Sf=None, n=10, n2=3, betas=None, freqs=None, cf=None):
    '''
    Initialize new filterbank
    topology: if specified (either as 'parallel' or 'series'), generates
      filterbank with filters taking on params and in the specified topology
    params: list of three lists of Ap's, bp's, Bu's, all of same length

    ALSO HAVE IN TERMS OF F NOT JUST BETA
    '''
    self.filters = []

    if filters is not None:
      for fil in filters:
        self.filters += [fil]
    else:
      args = {'type':type,
            'ir':ir, 'tf':tf,
            'coeffs':coeffs, 'roots':roots,
            'Ap':Ap, 'bp':bp, 'Bu':Bu, 'gain_const':gain_const, 'peak_magndb':peak_magndb,
            'Bpeak':Bpeak, 'fpeak':fpeak, 'phiaccum':phiaccum, 'Nbeta':Nbeta, 'Nf':Nf, 'Qerb':Qerb, 'ERBbeta':ERBbeta, 'ERBf':ERBf, 'Qn':Qn, 'Qn2':Qn2, 'BWndBbeta':BWndBbeta, 'BWndBbeta':BWndBf, 'BWn2dBbeta':BWn2dBbeta, 'BWn2dBf':BWn2dBf, 'Sbeta':Sbeta, 'Sf':Sf, 'n':n, 'n2':n2, 'betas':betas, 'freqs':freqs, 'cf':cf}
      scalar_args = dict()
      vector_args = dict()
      for k, v in args.items():
        if v is not None:
          if np.ndim(v) == 0:
            scalar_args[k] = v
          else:
            vector_args[k] = np.array(v)

      if not vector_args:
        self.filters = [Filter(**scalar_args)]
      else:
        vec_lens = [len(v) for v in vector_args.values()]
        if len(set(vec_lens)) != 1:
          raise Exception('All inputs to FilterBank must be either arraylikes of the same length or scalars')
        if vec_lens[0] == 0:
          warnings.warn('Len 0 input array provided, meaning 0 filters created')
        for i in range(vec_lens[0]):
          self.filters += [Filter(**scalar_args, **{k:vector_args[k][i] for k in vector_args})]
    # self.filters is now fully defined

    # testing: test arrays of different lengths, including 0, 1, many with or without scalars (incl. len 0 and scalars, which should work but be rather silly)

    self.graph = RootedTree()
    l = len(self.filters)
    if l == 1:
      self.graph.add()
    elif l > 1:
      if topology is None:
        raise Exception('Parallel/series topology must be indicated when initializing with multiple filters')
      elif topology=='parallel':
        self.graph.add_parallel(l)
      elif topology=='series':
        self.graph.add_series(l)
      else:
        raise Exception('Invalid topology name')

    self._uid2graphid = {self.filters[i].uid:i+1 for i in range(len(self.filters))}
    if len(self.filters) != len(self._uid2graphid):
      raise Exception('Filter UIDs must in fact be unique')
    self._uid2graphid[-1] = 0 # technically breaks if UID gets set to -1...
    self._graphid2uid = {v:k for k, v in self._uid2graphid.items()}

  def __len__(self):
    return len(self.filters)

  def get_filter_from_uid(self, uid):
    '''
    Identify filter from uid
    uid - uid
    '''
    return self.filters[self._uid2graphid[uid]-1]

  def get_source_uid(self, uid):
    return self._graphid2uid[self.graph.parent[self._uid2graphid[uid]]]

  def get_uids_fed_into(self, uid):
    return [self._graphid2uid[i] for i in self.graph.child[self._uid2graphid[uid]]]

  def add(self, filter, source=None, source_uid=-1):
    '''
    Add filter 'filter' to network, with source provided by 'source' (or can be given by UID).
    Source must be in network already.
    Supports various topologies.
    source_uid = -1 implies feed from input
    '''
    self.filters += [filter]
    graph_id = len(self.filters)
    self._uid2graphid[filter.uid] = graph_id
    if source is not None:
      graph_parentid = self._uid2graphid[source.uid]
    elif source_uid is not None:
      graph_parentid = self._uid2graphid[source_uid]
    else:
      raise Exception('Invalid source')
    self.graph.add(source_id=graph_parentid)

  def process_signal(self, signal, method='default'):
    '''
    Take input signal, feed through all filters
    '''
    def process(filter_input, graph_id):
      fil = self.get_filter_from_uid(self._graphid2uid[graph_id])
      return fil.solve(filter_input, method=method)

    outputs = self.graph.propagate_down(signal, process)
    return OutputSignals(outputs, self.graph)

  def bode_plot(self, num_samples=1000, freqs=None, peak_magndb=1, custom_title='Bode plot', show=True):
    if freqs is None:
      freqs = np.linspace(0.1, max(fil.get_computed_chars()['Bpeak'] for fil in self.filters)+1, num_samples)
    fils = [fil.bode_plot(freqs=freqs, peak_magndb=peak_magndb, show=False)+[fil.uid] for fil in self.filters]

    if show:
      fig, (ax1, ax2) = plt.subplots(2, 1)
      fig.suptitle(custom_title)

      for fil in fils:
        xaxis, magn, phase, uid = fil
        ax1.semilogx(xaxis, magn, label=uid) # magn in db
        ax2.semilogx(xaxis, phase, label=uid) # phase in cycles
      ax1.xaxis.set_major_locator(locator=matplotlib.ticker.LogLocator(subs=(1, 2, 5)))
      ax1.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f')) # pick better formatter
      ax1.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
      ax1.set_ylabel('Magnitude (dB)')
      # ax1.legend()
      ax2.xaxis.set_major_locator(locator=matplotlib.ticker.LogLocator(subs=(1, 2, 5)))
      ax2.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f')) # pick better formatter
      ax2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
      ax2.set_ylabel('Phase (radians)')
      ax2.set_xlabel('Normalized frequency')
      # ax2.legend()
      plt.show()
    return fils

if __name__ == "__main__":
  isignal = Signal.linear_chirp(f_init=1, f_final=10, fs=1000, num_samples=10000)
  # isignal = Signal(data=[1])
  isignal.plot()

  fb = FilterBank(topology='series', Ap=0.2, bp=[0.5, 1, 1.5], Bu=[3, 4, 5])
  fb.bode_plot(samples=np.linspace(0.05, 2.5, 10000))

  u = fb.filters[0].uid
  print(fb.get_filter_from_uid(u).uid == u) # True

  os = fb.process_signal(isignal, method='tf')
  f1 = fb.filters[0]

  # isignal = Signal.from_function(mode='t', func=(lambda t: np.sin(60*t)), fs=200, n=1000)
  # isignal *= [(1100-i)*(2+np.sin(i/35))/1000 for i in range(1000)]
  # isignal.plot()

  # f1 = Filter(Bpeak=1.5, Ncyc=11.1, phiaccum=3.5)
  # f1.bode_plot()
  # s1 = f1.solve(isignal, method='tf')
  # s1.plot()
  # for s in os.outsignals:
  #   s.plot()




  #########

  # isig = Signal.from_function(func=np.sin, fs=100, n=1000)
  # plt.plot(isig.timestamps, isig['t'])
  # # isig.plot()
  # num = 4
  # bplist = np.linspace(0.5, 1.5, num)
  # fs = FilterBank(topology='parallel', Ap=[0.1 for _ in range(4)], bp=bplist, Bu=[3 for _ in range(4)])
  # # fs.bode_plot()
  # osigs = fs.process_signal(isig)
  # # osigs.correlogram()
  # # f1, f2, f3 = fs.filters
  # for sig in osigs.outsignals:
  #   sig /= max(sig['t'])
  #   plt.plot(sig.timestamps, sig['t'])
  # #   sig.plot()
  # # s1, s2, s2 = osigs.signals
  # # s1 /= max(s1['t'])
  # # s2 /= max(s2['t'])
  # # plt.plot(s1.timestamps, s1['t'])
  # # plt.plot(s2.timestamps, s2['t'])
  # plt.show()

  # # isig = Signal.from_function(func=np.sin, fs=100, n=1000)
  # # plt.plot(isig.timestamps, isig['t'])
  # # # isig.plot()
  # # num = 4
  # # fs = FilterBank(topology='parallel', )
  # # # fs.bode_plot()
  # # osigs = fs.process_signal(isig)
  # # # osigs.correlogram()
  # # # f1, f2, f3 = fs.filters
  # # for sig in osigs.signals:
  # #   sig /= max(sig['t'])
  # #   plt.plot(sig.timestamps, sig['t'])
  # # #   sig.plot()
  # # # s1, s2, s2 = osigs.signals
  # # # s1 /= max(s1['t'])
  # # # s2 /= max(s2['t'])
  # # # plt.plot(s1.timestamps, s1['t'])
  # # # plt.plot(s2.timestamps, s2['t'])
  # # plt.show()

  # print(np.array([1, 2, 3]))

  # isig = Signal.from_function(func=np.sin, fs=100, n=1000)
  # plt.plot(isig.timestamps, isig['t'])
  # # isig.plot()
  # num = 4
  # bplist = np.linspace(0.5, 1.5, num)
  # fs = FilterBank(topology='parallel', Ap=[0.1 for _ in range(4)], bp=bplist, Bu=[3 for _ in range(4)])
  # # fs.bode_plot()
  # osigs = fs.process_signal(isig)
  # # osigs.correlogram()
  # # f1, f2, f3 = fs.filters
  # for sig in osigs.outsignals:
  #   sig /= max(sig['t'])
  #   plt.plot(sig.timestamps, sig['t'])
  # #   sig.plot()
  # # s1, s2, s2 = osigs.signals
  # # s1 /= max(s1['t'])
  # # s2 /= max(s2['t'])
  # # plt.plot(s1.timestamps, s1['t'])
  # # plt.plot(s2.timestamps, s2['t'])
  # plt.show()