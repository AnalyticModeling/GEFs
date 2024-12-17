import warnings
import numpy as np
import scipy.fft
import scipy.signal
import scipy.special
import scipy.integrate
import scipy.optimize
# import matplotlib.pyplot as plt

# from Signal import Signal
import helpers

class AbstractFilter():
  def __init__(self, uid, tf, ir, chars, betas):
    self.uid = uid
    self.tf = tf
    self.ir = ir
    self.betas = betas
    self.chars = chars

  def solve_fourier_tf(self, lendata, fftdata, freqs):
    X = fftdata
    Y = self.tf(1j*freqs)*X
    y = scipy.fft.irfft(Y, lendata)
    return y

  def solve_convolve_ir(self, data, timestamps):
    x = data
    ir = [self.ir(t) for t in timestamps]
    y = scipy.signal.convolve(ir, x)[:len(data)]
    return y

# double-check
  def _approx_func_ir_from_tf(self, tf, t, norm='backward', num_samples=10000, argmax_f=1):
    num_samples = -((-num_samples-1)//2)
    fs, freq_step = np.linspace(0, 1.5*argmax_f, num_samples, retstep=True)
    y = tf(fs)
    tot = y[0].real
    for idx in range(1, num_samples):
      tot += 2 * (y[idx] * np.exp(2j*np.pi*idx*t*freq_step)).real

    if norm=='backward':
      return tot/num_samples
    elif norm=='ortho':
      return tot/num_samples**0.5
    elif norm=='forward':
      return tot
    else:
      raise Exception('get the scipy error for incorrect normalization mode')

# double-check
  def _approx_func_tf_from_ir(self, ir, f, norm='backward', num_samples=10000, fs=100, time_span=None):
    if fs is not None:
      if time_span is not None:
        raise Exception('only provide num_samples and one of fs or time_span')
      time_span = (num_samples-1)/fs
    elif time_span is not None:
      fs = (num_samples-1)/time_span
    else:
      raise Exception('no input :(')

    x = [ir(i/fs) for i in range(num_samples)]
    tot = 0
    for idx in range(num_samples):
      tot += (x[idx] * np.exp(-2j*np.pi*idx*f/fs)).real

    if norm=='backward':
      return tot*num_samples
    elif norm=='ortho':
      return tot*num_samples**0.5
    elif norm=='forward':
      return tot
    else:
      raise Exception('get the scipy error for invalid normalization mode')

class Arbitrary(AbstractFilter):
  def __init__(self, uid, tf=None, ir=None, betas=None):
    if betas is None:
      betas = np.geomspace(0.01, 10, 10000)
    self.init_with_tf = (tf is not None)
    chars = helpers.computedfiltercharacteristics(tfunc=tf, betas=betas)

    if tf is not None:
      if ir is not None:
        raise Exception('Only one of tf and ir should be used to initialize Filter')
      def approx_ir(t, norm='backward', num_samples=10000, argmax_f=1):
        return self._approx_func_ir_from_tf(tf, t, norm=norm, num_samples=num_samples, argmax_f=chars['Bpeak'])
      super().__init__(uid, tf, approx_ir, chars, betas)
    elif ir is not None:
      def approx_tf(f, norm='backward', num_samples=10000, fs=100, time_span=None):
        return self._approx_func_tf_from_ir(ir, f, norm=norm, num_samples=num_samples, fs=fs, time_span=time_span)
      # auto compute length of time interval?
      super().__init__(uid, approx_tf, ir, chars, betas)
    else:
      raise Exception('At least one of tf and ir should be used to initialize Filter')

class Rational(AbstractFilter):
  def __init__(self, uid, coeffs=None, roots=None, betas=None):
    if betas is None: betas = np.geomspace(0.01, 10, 10000)
    if coeffs is not None:
      self.coeffs = coeffs
      self.order, self.roots, tf = self._given_coeffs(*coeffs)
    elif roots is not None:
      self.roots = roots
      self.order, self.coeffs, tf = self._given_roots(*roots)
    else:
      raise Exception('Invalid rational transfer function filter initialization')

    chars = helpers.computedfiltercharacteristics(tfunc=tf, betas=betas)
    def approx_ir(t, norm='backward', num_samples=10000, argmax_f=1):
        return self._approx_func_ir_from_tf(tf, t, norm=norm, num_samples=num_samples, argmax_f=chars['Bpeak'])
    super().__init__(uid, tf, approx_ir, chars, betas)

  def _given_coeffs(self, numer, denom):
    order = max(len(numer), len(denom))-1
    N = np.polynomial.polynomial.Polynomial(numer)
    D = np.polynomial.polynomial.Polynomial(denom)
    roots = [N.roots().tolist(), D.roots().tolist()]
    tf = lambda s: N(s)/D(s)
    return [order, roots, tf]

  def _given_roots(self, zeros, poles):
    order = max(len(zeros), len(poles))
    coeffs = [np.polynomial.polynomial.polyfromroots(zeros).tolist(), np.polynomial.polynomial.polyfromroots(poles).tolist()]
    N = np.polynomial.polynomial.Polynomial(coeffs[0])
    D = np.polynomial.polynomial.Polynomial(coeffs[1])
    tf = lambda s: N(s)/D(s)
    return [order, coeffs, tf]

  def solve_diffeq(self, data, fs=1):
    x = data

    x_coeffs, y_coeffs = self.coeffs
    timestamps = [i/fs for i in range(len(data))]
    x_derivatives = [x]

    for _ in range(len(x_coeffs)-1):
      x_derivatives += [np.gradient(x_derivatives[-1], timestamps, edge_order=2)]
    def nth_d_of_x(t, n):
      return np.interp(t, [j/fs for j in range(len(x))], x_derivatives[n])

    def update(t, y):
      LHS = sum(x_coeffs[i]*nth_d_of_x(t, i) for i in range(len(x_coeffs)))
      RHS = sum(y[i]*y_coeffs[i] for i in range(len(y)))
      y_last = (LHS - RHS) / y_coeffs[-1]
      dydt = [y[i+1] for i in range(len(y)-1)] + [y_last]
      return dydt
    y0 = [0 for _ in range(len(y_coeffs)-1)]
    return scipy.integrate.solve_ivp(update, (0, timestamps[-1]), y0, t_eval=timestamps).y[0]

class Parameterized(AbstractFilter):
  def __init__(self, uid, type='P', Ap=None, bp=None, Bu=None, gain_const=None, peak_magndb=None, Bpeak=None, phiaccum=None, Nbeta=None, Qerb=None, ERBbeta=None, Qn=None, Qn2=None, BWndBbeta=None, BWn2dBbeta=None, Sbeta=None, n=10, n2=3, betas=None):
    self.type = type
    if gain_const is None and peak_magndb is None: peak_magndb = 0
    self.gain_const = gain_const
    self.peak_magndb = peak_magndb
    if betas is None: betas = np.geomspace(0.01, 10, 10000)

    has_params = any(param is not None for param in [Ap, bp, Bu])
    chars = dict()
    if any(characteristic is not None for characteristic in [Bpeak, phiaccum, Nbeta, Qerb, Qn, Qn2, Sbeta]):
      chars = {'Bpeak':Bpeak, 'phiaccum':phiaccum, 'Nbeta':Nbeta, 'Qerb':Qerb, 'ERBbeta':ERBbeta, 'Qn':Qn, 'Qn2':Qn2, 'BWndBbeta':BWndBbeta, 'BWn2dBbeta':BWn2dBbeta, 'Sbeta':Sbeta}
      for k in ['Bpeak', 'phiaccum', 'Nbeta', 'Qerb', 'ERBbeta', 'Qn', 'Qn2', 'BWndBbeta', 'BWn2dBbeta', 'Sbeta']:
        if chars[k] == None:
          del chars[k]
    if has_params and any(param is None for param in [Ap, Bu]):
      raise Exception('Model parameters Ap and Bu should be provided')
    if chars and sum(c in chars for c in ['phiaccum','Nbeta','Qerb','Qn','Qn2','Sbeta']) != 2:
      raise Exception('Invalid set of filter characteristics provided')

    if has_params: # add more checks here
      if bp is None: bp = 1
      if Ap < 0: raise Exception('Ap should not be negative')
      if bp < 0: raise Exception('bp should not be negative')
      self.params = {'Ap':Ap, 'bp':bp, 'Bu':Bu}
      tf, chars = self._given_p(Ap, bp, Bu, betas)
      self.orig_chars = chars
    else:
      for k in ['Bpeak', 'Qerb', 'ERBbeta', 'Qn', 'Qn2', 'BWndBbeta', 'BWn2dBbeta', 'Sbeta']:
        if k in chars and chars[k] < 0:
          raise Exception(f'{k} should not be negative')
      self.orig_chars = {k:v for k, v in chars.items()}
      if self.orig_chars is None: raise Exception('shouldbe impossible to get here')
      tf, self.params, chars = self._given_c(chars, n, n2, betas)

    def double_fact(x):
      if abs(x-round(x)) < 1e-10:
        return 2**(x-1)*scipy.special.gamma(x)
      else:
        return (2/np.pi)**0.5*2**(x-1)*scipy.special.gamma(x)

    def param_ir(t):
      # warnings.warn('Parametrized approximate impulse response may produce incorrect values for very small t')
      Ap = self.params['Ap']
      bp = self.params['bp']
      Bu = self.params['Bu']
      return (1/(double_fact(Bu)*bp**(Bu-0.5))) * np.exp(-Ap*t) * (t**(Bu-1/2)) * scipy.special.jn(Bu-1/2, t*bp)
      # write cosine approximation somewhere

    super().__init__(uid, tf, param_ir, chars, betas)

  def _given_p(self, Ap, bp, Bu, betas):
    p = -Ap + 1j*bp
    if self.type == 'P':
      temp_tf = lambda s: ((s - p) * (s - p.conjugate()))**(-Bu)
    elif self.type == 'V':
      temp_tf = (lambda s: ((s + Ap) / ((s - p) * (s - p.conjugate()))**(Bu+1)))

    else:
      raise Exception('Filter type not defined')

    if self.gain_const is not None:
      C = self.gain_const
    else:
      if self.peak_magndb is None:
        raise Exception('Should be impossible to get here')
      curr_peak_magn = max(abs(np.array([temp_tf(x*1j) for x in betas])))
      peak_magn = np.power(10, self.peak_magndb/20)
      C = peak_magn/curr_peak_magn
    tf = lambda s: C * temp_tf(s)

    chars = helpers.computedfiltercharacteristics(tfunc=tf, betas=betas)
    return [tf, chars]

  def _given_c(self, chars, n, n2, betas):
    # get closed forms?
    params = helpers.chars2params(chars, n=n, n2=n2)
    tf, cs = self._given_p(params['Ap'], params['bp'], params['Bu'], betas)
    return [tf, params, cs]

  def sharp_approximation(self, n=10, n2=3):
    Ap = self.params['Ap']
    if Ap > 0.3:
      warnings.warn('Sharp filter approximation may not be valid for Ap > 0.3')
    bp = self.params['bp']
    Bu = self.params['Bu']
    p = -Ap + 1j*bp
    tf = lambda s: (s - p)**(-Bu)
    return [tf, helpers.sharpfiltercharacteristics(Ap, bp, Bu, n=n, n2=n2)]

  def solve_as_rational_tf(self, data, fs=1):
    x = data

    if isinstance(self.params['Bu'], (int, np.integer)):
      Bu = self.params['Bu']
    else:
      warnings.warn('Rounding Bu to nearest integer')
      Bu = round(self.params['Bu'])

    y_coeffs = np.polynomial.polynomial.polypow([self.params['Ap']**2+self.params['bp']**2, 2*self.params['Ap'], 1], Bu)
    timestamps = [i/fs for i in range(len(data))]

    def update(t, y):
      LHS = np.interp(t, timestamps, x)
      RHS = sum(y[i]*y_coeffs[i] for i in range(len(y)))
      y_last = (LHS - RHS) # /y_coeffs[-1] is just /1
      dydt = [y[i+1] for i in range(len(y)-1)] + [y_last]
      return dydt

    y0 = [0 for _ in range(len(y_coeffs)-1)]
    Ptypesol = scipy.integrate.solve_ivp(update, (0, timestamps[-1]), y0, method='RK45', t_eval=timestamps).y[0]

    if self.type == 'V':
      output = np.gradient(Ptypesol, timestamps) + self.params['Ap']*Ptypesol
    else:
      output = Ptypesol
    return output

  def solve_fractional_diffeq(self, data, fs=1):
    Ap = self.params['Ap']
    bp = self.params['bp']
    p = -Ap+1j*bp
    Bu = self.params['Bu']

    delta_t = 1/fs
    N = len(data)

    enegpTuT = [np.exp(-p*Tscaled*delta_t) * data[Tscaled] for Tscaled in range(N)]
    inner_int = []
    for tauscaled in range(N):
      row = [((tauscaled-Tscaled)*delta_t)**(Bu-1) * enegpTuT[Tscaled] for Tscaled in range(tauscaled+1)]
      inner_int += [sum(row)*delta_t]
    e2ibptau_etc = [np.exp(2j*bp*tauscaled*delta_t) * inner_int[tauscaled] for tauscaled in range(N)]
    outer_int = []
    for tscaled in range(N):
      row = [((tscaled-tauscaled)*delta_t)**(Bu-1) * e2ibptau_etc[tauscaled] for tauscaled in range(tscaled+1)]
      outer_int += [sum(row)*delta_t]
    Ptypesol_approx = [(1/scipy.special.gamma(Bu))**2 * np.exp(p.conjugate()*(tscaled*delta_t)) * outer_int[tscaled] for tscaled in range(N)]
    Ptypesol = [z.real for z in Ptypesol_approx]

    if self.type == 'V':
      deriv = np.gradient(Ptypesol, [i/fs for i in range(len(data))])
      output = np.array([deriv[i]+self.params['Ap']*Ptypesol[i] for i in range(len(data))])
    else:
      output = Ptypesol
    return output