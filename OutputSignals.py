import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from Signal import Signal

class OutputSignals:
  '''
  Output of processing Signal through a FilterBank.
  '''
  def __init__(self, all_signals, graph):
    '''
    Users should not be personally instantiating these objects.
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
    self._uid2graphid[-1] = 0 # technically breaks if UID gets set to -1... (which is why this should be treated as protected)
    self._graphid2uid = {v:k for k, v in self._uid2graphid.items()}

  def __len__(self):
    return len(self.outsignals)

  def get_signal_from_uid(self, uid):
    '''
    Get Signal in OutputSignals given its UID

    Arguments:
      uid: UID of Signal to be identified
    '''
    return self.outsignals[self._uid2graphid[uid]-1]

  def get_source_uid(self, uid):
    '''
    Given a UID corresponding to a Signal under consideration, provides the UID \
      of the original Signal that was run through a Filter to generate the Signal under consideration.

    Returning -1 means the Signal fed into the Filter was the original `self.insignal`.

    Arguments:
      uid: UID of Signal whose source is desired
    '''
    return self._graphid2uid[self.graph.parent[self._uid2graphid[uid]]]

  def correlogram(self, custom_title='Correlogram'):
    '''
    Draws correlogram, i.e. the plot of all possible pairs of correlations between Signals in OutputSignals object.

    Arguments:
      custom_title: Optional title of plot. Default is 'Correlogram'.
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

    for i in range(n-1):
      axs[-1][i].set_xlabel(f'Signal {self.outsignals[i].uid}')
      axs[i][0].set_ylabel(f'Signal {self.outsignals[i+1].uid}')

    fig.subplots_adjust(wspace=0, hspace=0)
    fig.suptitle(custom_title)
    fig.tight_layout()

    plt.show()

  def autocorrelates(self, custom_title='All autocorrelates'):
    '''
    Simultaneously plots all autocorrelations of Signals in OutputSignals object.

    Arguments:
      custom_title: Optional title of plot. Default is 'All autocorrelates'.
    '''
    n = len(self.outsignals)
    fig, axs = plt.subplots(n)
    for i in range(n):
      subgraph = axs[i]
      subgraph.plot(range(self.signal_length), self.outsignals[i].autocorrelate())
      subgraph.set_ylabel(f'Signal {self.outsignals[i].uid}')

    fig.subplots_adjust(wspace=0, hspace=0)
    fig.suptitle(custom_title)

    plt.show()

  def correlate_with(self, signal, custom_title='Correlations'):
    '''
    Simultaneously plots all correlations of the Signal `signal` with \
      all Signals in this OutputSignals object.

    Arguments:
      custom_title: Optional title of plot. Default is 'Correlations'.
    '''
    n = len(self.outsignals)
    fig, axs = plt.subplots(n)
    for i in range(n):
      subgraph = axs[i]
      subgraph.plot(range(-self.signal_length+1, self.signal_length), signal @ self.outsignals[i])
      subgraph.set_ylabel(f'Signal {self.outsignals[i].uid}')

    fig.subplots_adjust(wspace=0, hspace=0)
    fig.suptitle(custom_title)

    plt.show()