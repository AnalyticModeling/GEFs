import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from Signal import Signal

class OutputSignals:
  def __init__(self, all_signals, graph):
    '''
    Users probably shouldn't make their own instances
    '''
    if not all(isinstance(s, Signal) for s in all_signals): raise TypeError('All inputs to OutputSignals should be Signal objects')
    self.insignal = all_signals[0]
    self.signal_length = len(self.insignal)
    self.signal_fs = self.insignal.fs
    self.outsignals = all_signals[1:]
    self.graph = graph
    self._uid2graphid = {self.outsignals[i].uid:i+1 for i in range(len(self.outsignals))}
    if len(self.outsignals) != len(self._uid2graphid):
      raise Exception('Signal UIDs must in fact be unique')
    self._uid2graphid[-1] = 0 # technically breaks if UID gets set to -1...
    self._graphid2uid = {v:k for k, v in self._uid2graphid.items()}

  def __len__(self):
    return len(self.outsignals)

  def get_signal_from_uid(self, uid):
    return self.outsignals[self._uid2graphid[uid]-1]

  def get_source_uid(self, uid):
    return self._graphid2uid[self.graph.parent[self._uid2graphid[uid]]]

  def correlogram(self, custom_title='Correlogram'):
    '''
    draw correlogram
    Y AXIS SHOULD BE LABELED OR ELSE NO CLUE HOW HIGH IT IS
    '''
    n = len(self.outsignals)
    fullxaxis = [i/self.signal_fs for i in range(-self.signal_length+1, self.signal_length)]
    fig, axs = plt.subplots(n-1, n-1)

    for i in range(n-1):
      for j in range(n-1):
        if i < j:
          fig.delaxes(axs[i][j])
          continue

        fst = self.outsignals[i+1]
        snd = self.outsignals[j]
        corr = fst @ snd
        subgraph = axs[i][j]
        subgraph.plot(fullxaxis, corr)
        if i != n-1:
          subgraph.set_xticks([])
        # if j != 0:
        #   subgraph.set_yticks([])
        # axs[i][j].set_axis_off()

    for i in range(n-1):
      axs[-1][i].set_xlabel(f'Signal {self.outsignals[i].uid}')
      axs[i][0].set_ylabel(f'Signal {self.outsignals[i+1].uid}')

    fig.subplots_adjust(wspace=0, hspace=0)
    fig.suptitle(custom_title)
    fig.tight_layout()

    plt.show()

  def autocorrelates(self, custom_title='Autocorrelates'):
    n = len(self.outsignals)
    fig, axs = plt.subplots(n)
    for i in range(n):
      subgraph = axs[i]
      # print('testing', len(self.outsignals[i].autocorrelate()))
      # print(len(range(self.signal_length)))
      subgraph.plot(range(self.signal_length), self.outsignals[i].autocorrelate())
      subgraph.set_ylabel(f'Signal {self.outsignals[i].uid}')

    fig.subplots_adjust(wspace=0, hspace=0)
    fig.suptitle(custom_title)

    plt.show()

  def correlate_with(self, signal, custom_title='Correlations'):
    n = len(self.outsignals)
    fig, axs = plt.subplots(n)
    for i in range(n):
      subgraph = axs[i]
      subgraph.plot(range(-self.signal_length+1, self.signal_length), signal @ self.outsignals[i])
      subgraph.set_ylabel(f'Signal {self.outsignals[i].uid}')

    fig.subplots_adjust(wspace=0, hspace=0)
    fig.suptitle(custom_title)

    plt.show()

if __name__ == "__main__":
  # isignal = Signal.linear_chirp(f_init=1, f_final=10, fs=1000, n=10000)
  # # d = isignal['t']
  # isignal.plot()
  # plt.plot([i/isignal.fs for i in range(len(d))], d)
  # plt.show()

  ########

  sr = 100
  num = 1000
  isignal = Signal.linear_chirp(f_init=1, f_final=1, fs=sr, n=num)
  s1 = Signal.linear_chirp(f_init=5, f_final=0.5, fs=sr, n=num)
  s1 *= (1.5+np.sin([i/sr for i in range(num)]))
  s2 = Signal.linear_chirp(f_init=2, f_final=1, fs=sr, n=num)
  s2 *= (1.5+np.sin([i/sr+2 for i in range(num)]))
  s3 = Signal.linear_chirp(f_init=5, f_final=2, fs=sr, n=num)
  s3 *= (1.5+np.sin([i/sr+1 for i in range(num)]))
  # xaxis = [i/sr for i in range(num)]
  # plt.plot(xaxis, s1['t'])
  # plt.plot(xaxis, s2['t'])
  # plt.plot(xaxis, s3['t'])
  # # u1, l1 = s1.envelope_analytic()
  # # plt.plot(xaxis, u1)
  # # plt.plot(xaxis, l1)
  # # u2, l2 = s2.envelope_analytic()
  # # plt.plot(xaxis, u2)
  # # plt.plot(xaxis, l2)
  # plt.show()

  # print(isignal.uid, s1.uid, s2.uid)

  # fullxaxis = [i/sr for i in range(-num+1, num)]
  # plt.plot(fullxaxis, s1 @ s2)
  # plt.show()

  os = OutputSignals([isignal, s1, s2, s3], None)
  # os.add(s1)
  # os.add(s2)
  # os.add(s3)
  os.correlogram()
  os.autocorrelates()
  # print(len(isignal))
  # print(len(s1))
  os.correlate_with(isignal)



