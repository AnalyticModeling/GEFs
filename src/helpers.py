import numpy as np
# import scipy as sp
import scipy.optimize
import scipy.special
import matplotlib.pyplot as plt

def mag2db(x):
  return 20*np.log10(x)

def match_lengths(*args, exception_msg='Input list lengths do not match'):
  arg_copy = []
  list_exist = False
  list_len = -1
  for arg in args:
    if np.ndim(arg) != 0:
      list_exist = True
      list_len = len(arg)
      break

  if list_exist:
    for arg in args:
      if np.ndim(arg) == 0:
        arg_copy += [[arg for _ in range(list_len)]]
      elif len(arg) != list_len:
        raise Exception(exception_msg)
      else:
        arg_copy += [arg]
    return arg_copy
  else:
    return [[arg] for arg in args] # returns a list of len 1 if all inputs are scalars

def chars2consts(Psi, n=10, n2=3):
  '''
  The inputs and outputs of this function are always in terms of normalized frequency.
  '''
  a = 0.418
  b = 1.02

  if 'Bpeak' in Psi:
    Bpeak = Psi['Bpeak']
  else:
    Bpeak = 1
  bp = Bpeak

  if ('ERBbeta' in Psi): Psi['Qerb'] = Bpeak/Psi['ERBbeta']
  if ('BWndBbeta' in Psi): Psi['Qn'] = Bpeak/Psi['BWndBbeta']
  if ('BWn2dBbeta' in Psi): Psi['Qn2'] = Bpeak/Psi['BWn2dBbeta']

  if ('Qn2' in Psi):
    if not (('Qn' in Psi) and ('Qn2' in Psi)): raise Exception
    Qn = Psi['Qn']
    Qn2 = Psi['Qn2']
    Bu = scipy.optimize.fsolve(lambda x: -Qn/Qn2 + ((10 ** (n2 / 10 / x) - 1)/(10 ** (n / 10 / x) - 1)) ** -0.5, 3)
    Ap = bp / 2 / Qn / (10 ** (n / 10 / Bu) - 1) ** 0.5
  elif 'phiaccum' in Psi:
    phiaccum = Psi['phiaccum']
    Bu = 2 * phiaccum
    if 'Nbeta' in Psi: Ap = phiaccum / np.pi / Psi['Nbeta']
    if 'Qerb' in Psi: Ap = (Bpeak / (np.pi)**0.5 / Psi['Qerb']) * (scipy.special.gamma(Bu) / scipy.special.gamma(Bu - 0.5))
    if 'Qn' in Psi: Ap = (Bpeak / 2 / Psi['Qn']) * ((10 ** (n / 20 / phiaccum))-1)
    if 'Sbeta' in Psi: Ap = ((40 * np.pi / np.log(10)) * (phiaccum / Psi['Sbeta'])) ** 0.5
  elif 'Nbeta' in Psi:
    # print('N in Psi')
    Nbeta = Psi['Nbeta']
    if 'Qerb' in Psi:
      Bu = np.exp(b/a) * (Bpeak * Nbeta / Psi['Qerb']) ** (1/a)
    elif 'Qn' in Psi:
      Bu = scipy.optimize.fsolve(lambda x: ((-1 + 10 ** (n / 10 / x)) * np.pi / x) - (Psi['Qn'] / Bpeak / Nbeta), 3)
    elif 'Sbeta' in Psi:
      Bu = 80 * np.pi**2 * Nbeta**2 / np.log(10) / Psi['Sbeta']
    else: raise Exception('Should be impossible to have no value for Bu')
    Ap = Bu / 2 / np.pi / Nbeta
  else:
    if 'Qn' in Psi and 'Qerb' in Psi:
      Qn, Qerb = Psi['Qn'], Psi['Qerb']
      Bu = scipy.optimize.fsolve(lambda x: -Qerb/Qn + (2 / np.pi**0.5) * (scipy.special.gamma(x) / scipy.special.gamma(x - 0.5)) * (-1 + 10 ** (n / 10 / x))**0.5, 3)
      Ap = bp / 2 / Qn / (-1 + 10 ** (n / 10 / Bu))**0.5
    elif 'Qn' in Psi and 'Sbeta' in Psi:
      Qn, Sbeta = Psi['Qn'], Psi['Sbeta']
      Bu = scipy.optimize.fsolve(lambda x: -27.942143326758714 + 34.74355855226014 * x * (-1 + 10 ** (n / 10 / x)) / bp**2, 3)
      Bu = scipy.optimize.fsolve(lambda x: -Sbeta/Qn**2 + (80 / np.log(10)) * x * (-1 + 10 ** (n / 10 / x)) / bp**2, 3)
      Ap = (20 * Bu / np.log(10) / Sbeta)**0.5
    elif 'Qerb' in Psi and 'Sbeta' in Psi:
      Qerb, Sbeta = Psi['Qerb'], Psi['Sbeta']
      Bu = scipy.optimize.fsolve(lambda x: -Sbeta/Qerb**2 + (20 * np.pi / np.log(10)) * scipy.special.gamma(x - 0.5)**2 / scipy.special.gamma(x)**2 / bp**2, 3)
      Ap = (20 * Bu / np.log(10) / Sbeta)**0.5
    else: raise Exception('Insufficient number of filter characteristics')

  return {'Ap': float(Ap), 'bp': float(bp), 'Bu': float(Bu)}

def computedfiltercharacteristics(tfunc=None, data=None, betas=None, n=10., n2=3.):
  '''
  Given filter characteristics, computes approximate parameters from transfer function
  tfunc: transfer function
  sample: points to sample along
  n: bandwidth parameter
  '''
  def find_nearest_index(arr, x):
    return np.argmin(abs(arr-x))

  # NOTE_TO_SELF: calculations done in cycles (so period = 1)

  if betas is None: betas = np.geomspace(0.01, 10, 10000)
  if data is None:
    if tfunc is None: raise Exception('Transfer function should be provided in either function or list form')
    data = np.array([tfunc(x*1j) for x in betas])
  phases_raw = np.unwrap(np.angle(data))
  phases_cyc = phases_raw/(2*np.pi)

  magns_raw = abs(data)
  magns_db = 20*np.log10(magns_raw)
  Bpeak_index = np.argmax(magns_raw)
  Bpeak = betas[Bpeak_index]

  N = -min(np.gradient(phases_cyc, betas))
  phiaccum = (max(phases_cyc)-min(phases_cyc)) # is this phiaccum?
  def get_BWndBbeta(num):
    # maybe do this smarter?
    if Bpeak_index == 0:
      BWndBbeta_index = betas[find_nearest_index(magns_db, magns_db[Bpeak_index]-num)]-betas[0]
    elif Bpeak_index == len(betas)-1:
      BWndBbeta_index = betas[-1]-betas[find_nearest_index(magns_db, magns_db[Bpeak_index]-num)]
    elif magns_db[Bpeak_index] == float('inf'):
      raise Exception('Transfer function should not evaluate to an infinite magnitude')
    else:
      BWndBbeta_index = betas[Bpeak_index+1+find_nearest_index(magns_db[Bpeak_index+1:], magns_db[Bpeak_index]-num)] \
                        - betas[find_nearest_index(magns_db[:Bpeak_index], magns_db[Bpeak_index]-num)]
    return BWndBbeta_index
  BWndBbeta = get_BWndBbeta(n)
  Qn = Bpeak/BWndBbeta
  BWn2dBbeta = get_BWndBbeta(n2)
  Qn2 = Bpeak/BWn2dBbeta

  ERBbeta = np.trapz((magns_raw/magns_raw[Bpeak_index])**2, x=betas)
  Qerb = Bpeak/ERBbeta

  Sbeta = -np.gradient(np.gradient(magns_db, betas), betas)[Bpeak_index]

  # Qerb2N = Qerb/N
  # Qn2N = Qn/N
  # Qerb2Qn = Qerb/Qn

  return {'Bpeak':Bpeak,
          'Nbeta':N,
          'phiaccum':phiaccum,
          'Qerb':Qerb,
          'ERBbeta':ERBbeta,
          'Qn':Qn,
          'n':n,
          'BWndBbeta':BWndBbeta,
          'Qn2':Qn2,
          'n2':n2,
          'BWn2dBbeta':BWn2dBbeta,
          'Sbeta':Sbeta}

def sharpfiltercharacteristics(Ap, bp, Bu, n=10, n2=3):
  '''
  Given filter constants (Ap, bp, Bu, n, n2), computes filter characteristics with sharp-filter approximation.
  Assumes approximation holds (Filter.py checks if this is actually so)

  Attributes:
    Ap, bp, Bu: model parameters
    n, n2: bandwidth parameters (two will always be computed; the second can be ignored in the output if undesired)
  '''
  Bpeak = bp
  N = Bu/(2*np.pi*Ap)
  phiaccum = Bu/2
  BWndBbeta = 2*Ap * (-1 + 10**(n/(10*Bu)))**0.5
  Qn = bp/BWndBbeta
  BWn2dBbeta = 2*Ap * (-1 + 10**(n2/(10*Bu)))**0.5
  Qn2 = bp/BWn2dBbeta
  ERBbeta = (np.pi**0.5)*Ap*scipy.special.gamma(Bu-1/2)/scipy.special.gamma(Bu)
  Qerb = bp/ERBbeta
  Qerb2N = Qerb/N
  Qn2N = Qn/N
  Qerb2Qn = Qerb/Qn
  Sbeta = 20/np.log(10)*Bu/Ap**2
  return {'Bpeak':Bpeak,
          'Nbeta':N,
          'phiaccum':phiaccum,
          'Qerb':Qerb,
          'ERBbeta':ERBbeta,
          'Qn':Qn,
          'Qn2':Qn2,
          'BWndBbeta':BWndBbeta,
          'BWn2dBbeta':BWn2dBbeta,
          'Sbeta':Sbeta,
          'Qerb2N':Qerb2N,
          'Qn2N':Qn2N,
          'Qerb2Qn':Qerb2Qn}

def plot_with_arrow(x, y, xlabel=None, ylabel=None, custom_title=None):
  plt.plot(x, y)
  if custom_title is not None:
    plt.title(custom_title)
  if xlabel is not None:
    plt.xlabel(xlabel)
  if ylabel is not None:
    plt.ylabel(ylabel)
  arrow_length = 0.05
  dx = x[-1]-x[-2]
  dy = y[-1]-y[-2]
  theta = np.arctan2(dy, dx)
  plt.arrow(x[-2], y[-2], np.cos(theta)*arrow_length, np.sin(theta)*arrow_length, head_width=0.01, head_length=0.01, length_includes_head=True)
  plt.show()

def plot_2x1(xaxis, upper, lower, xlabel=None, upper_ylabel=None, lower_ylabel=None, custom_title=None):
  fig, (ax1, ax2) = plt.subplots(2, 1)

  ax1.plot(xaxis, upper)
  if upper_ylabel is not None:
    ax1.set_ylabel(upper_ylabel)

  ax2.plot(xaxis, lower)
  if lower_ylabel is not None:
    ax2.set_ylabel(lower_ylabel)

  if xlabel is not None:
    ax2.set_xlabel(xlabel)
  if custom_title is not None:
    fig.suptitle(custom_title)

  plt.show()