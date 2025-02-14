from Filter import Filter
from Signal import Signal
from OutputSignals import OutputSignals
from RootedTree import RootedTree

import warnings
import numpy as np
# import scipy as sp
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker

class FilterBank:
  '''
  Object storing multiple Filters in a graph structure to modeling being able to feed a Signal into a Filter and having the output feed into another Filter.
  '''
  def __init__(self, topology=None, filters=None, type='P', ir=None, tf=None, coeffs=None, roots=None, Ap=None, bp=None, Bu=None, gain_const=None, peak_magndb=None, Bpeak=None, fpeak=None, phiaccum=None, Nbeta=None, Nf=None, Qerb=None, ERBbeta=None, ERBf=None, Qn=None, Qn2=None, BWndBbeta=None, BWndBf=None, BWn2dBbeta=None, BWn2dBf=None, Sbeta=None, Sf=None, n=10, n2=3, betas=None, freqs=None, cf=None):
    '''
    Initialize new filterbank. Most arguments are the same as for `Filter` object, \
      though if two arguments are vectors, they must be the same length. Scalars will \
      be broadcast to the same same as vectors. If everything is a scalar, it is all \
      interpreted as a vector of length 1.

    'parallel' means Filters all take input from the same source. 'series' means \
      Filters take

    Arguments:
      topology: If specified, initializes filterbank with filters taking on params \
        and in the specified topology. There are two options:
        - 'parallel': all Filters take input from the same source
        - 'series': each Filter takes input from the output of the Filter before it in sequence
      filters: If a list of Filters already exists, these Filters can just be placed into \
        the FilterBank, skipping having to initialize everything
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
    Get Filter in FilterBank given its UID

    Arguments:
      uid: UID of Filter to be identified
    '''
    return self.filters[self._uid2graphid[uid]-1]

  def get_source_uid(self, uid):
    '''
    Given a UID corresponding to a Filter, provides the UID of the source of that Filter. \
      Returning -1 means no other Filter feeds into said Filter.

    Arguments:
      uid: UID of Filter whose source is desired
    '''
    return self._graphid2uid[self.graph.parent[self._uid2graphid[uid]]]

  def get_uids_fed_into(self, uid):
    '''
    Given a UID corresponding to a Filter, provides the UIDs of the 'children' \
      of that Filter (i.e. the Filters that accept the output of that Filter as input).

    Arguments:
      uid: UID of Filter whose children are desired
    '''
    return [self._graphid2uid[i] for i in self.graph.child[self._uid2graphid[uid]]]

  def add(self, filter, source=None, source_uid=-1):
    '''
    Add Filter 'filter' to network.

    Either the source filter can be directly provided or the source can be referred to by UID.

    A source UID of -1 corresponds to the Filter directly being fed input not from any other Filter.

    Source must be in network already.

    Arguments:
      filter: Filter to be added
      source: Source Filter
      source_uide: UID of source Filter. Default is -1
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

  def process_signal(self, signal, method=None):
    '''
    Given a Signal, feeds the Signal through the FilterBank and returns \
      an OutputSignals object. These are guaranteed to share a topology.

    Arguments:
      signal: Signal to processed by FilterBank
      method: Method used to process. See `Filter.solve()` for more details
    '''
    def process(filter_input, graph_id):
      fil = self.get_filter_from_uid(self._graphid2uid[graph_id])
      return fil.solve(filter_input, method=method)

    outputs = self.graph.propagate_down(signal, process)
    return OutputSignals(outputs, self.graph)

  def bode_plot(self, freqs=None, custom_title='Bode plot', show=True):
    '''
    Generate simultaneous Bode plots of all Filters in FilterBank.

    Returns list of quadruples of the form [x-axis (frequency) data, magnitudes (dB), phases (cycles), filter UID].

    See `Filter.bode_plot` for more details on arguments.
    '''
    if freqs is None:
      freqs = np.geomspace(0.01, 3*max(fil.get_computed_chars()['Bpeak'] for fil in self.filters), 10000)
    fils = [fil.bode_plot(freqs=freqs, show=False)+[fil.uid] for fil in self.filters]

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

      ax2.xaxis.set_major_locator(locator=matplotlib.ticker.LogLocator(subs=(1, 2, 5)))
      ax2.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f')) # pick better formatter
      ax2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
      ax2.set_ylabel('Phase (radians)')
      ax2.set_xlabel('Normalized frequency')

      plt.show()
    return fils