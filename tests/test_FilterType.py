import pytest
import numpy as np
from src.GEFs_core.Filter import Filter
from src.GEFs_core.FilterType import AbstractFilter, Arbitrary, Rational, Parameterized
from src.GEFs_core.Signal import Signal
from src.GEFs_core.helpers import chars2consts, sharpfiltercharacteristics

#======================
# Filters Fixture
#======================

@pytest.fixture(params=[
    # Filters parameters that will be used to create/initialize the filters. (can be extended to include more filters)
    {'type': 'Arbitrary', 'tf': lambda s: 1 / (1 + s + s**2)},
    {'type': 'Arbitrary', 'ir': lambda t: np.exp(-t)},
    {'type': 'Rational', 'coeffs': [[1], [1, 0.5, 0.25]]},
    {'type': 'Rational', 'roots': ([-1], [0.5 + 0.5j, 0.5 - 0.5j])},
    {'type': 'Parameterized', 'Ap': 0.1, 'bp': 1, 'Bu': 3, 'cf': 10.0},
    {'type': 'Parameterized', 'phiaccum': np.pi / 4, 'Qerb': 1.0},
])
def filter_under_test(request):
    """ Create/initialize various Filter objects for testing by requesting the parameters from the filters fixture."""
    params = request.param
    filter_type = params.pop('type')

    if filter_type == 'Arbitrary':
        if 'tf' in params:
            return Arbitrary(uid=1, tf=params['tf'])
        elif 'ir' in params:
            return Arbitrary(uid=2, ir=params['ir'])
        else:
            raise ValueError("Arbitrary filter must have either 'tf' or 'ir' defined.")
    elif filter_type == 'Rational':
        if 'coeffs' in params:
            return Rational(uid=3, coeffs=params['coeffs'])
        elif 'roots' in params:
            return Rational(uid=4, roots=params['roots'])
        else:
            raise ValueError("Rational filter must have either 'coeffs' or 'roots' defined.")
    elif filter_type == 'Parameterized':
        if 'Ap' in params and 'bp' in params and 'Bu' in params:
            return Parameterized(uid=5, Ap=params['Ap'], bp=params['bp'], Bu=params['Bu'])
        elif 'phiaccum' in params and 'Qerb' in params:
            return Parameterized(uid=6, phiaccum=params['phiaccum'], Qerb=params['Qerb'])
        else:
            raise ValueError("Parameterized filter must have either 'Ap', 'bp', and 'Bu' or 'Bpeak' and 'phiaccum' defined.")
    else:
        raise ValueError(f"Unknown filter type: {filter_type}")

#======================
# Signals Fixture
#======================

@pytest.fixture(params=[
    # Signals parameters that will be used to create/initialize the signals. (can be extended to include more signals)
    {
        'mode': 't',
        'func': lambda x: (1 + x / 100) * np.sin(x / 10),
        'num_samples': 100,
        'fs': 10
    },
    {
        'mode': 'f',
        'func': lambda x: x / 10,
        'num_samples': 20,
        'evenlen': False,
        'fs': 10
    }
])
def signal_under_test(request):
    """ Create/initialize various signals objects for testing by requesting the parameters from the signals fixture."""
    params = request.param
    mode = params.pop('mode')
    func = params.pop('func')
    num_samples = params.pop('num_samples')
    fs = params['fs']

    if mode == 't':
        return Signal.from_function(mode=mode, func=func, num_samples=num_samples, fs=fs)
    elif mode == 'f':
        return Signal.from_function(mode=mode, func=func, num_samples=num_samples, fs=fs)
    else:
        raise ValueError(f"Unknown signal mode: {mode}")

#======================
# Core Tests
#======================

# Test 1 - Passed
def test_abstract_filter_initialization_properties(filter_under_test):
    """Test initialization properties of AbstractFilter."""
    assert isinstance(filter_under_test.uid, int)
    assert callable(filter_under_test.tf)
    assert callable(filter_under_test.ir)
    assert isinstance(filter_under_test.chars, dict)
    assert isinstance(filter_under_test.betas, np.ndarray)

# Test 2 - Passed (maybe uid should be specified to support only int type? also it solves given ir and tf together but i don't know if it needs to be specified about this in the abstract class? )
def test_abstract_filter_init():
    uid = "filter1"
    tf = lambda freq: freq ** 2
    ir = lambda t: t + 1
    chars = {"Qerb": 1.0, "fake": "fake"}  # Only one characteristic and 1 invalid, but ok for abstract class
    betas = [0.5, 1.0]
    abstract_filter = AbstractFilter(uid, tf, ir, chars, betas)

    print(f"Abstract Filter UID: {abstract_filter.uid}")
    print(f"Transfer Function (tf(2)): {abstract_filter.tf(2)}")
    print(f"Impulse Response (ir(3)): {abstract_filter.ir(3)}")
    print(f"Characteristics (chars): {abstract_filter.chars}")
    print(f"Betas: {abstract_filter.betas}")

    assert abstract_filter.uid == uid
    assert abstract_filter.tf(2) == 4
    assert abstract_filter.ir(3) == 4
    assert abstract_filter.chars == chars
    assert abstract_filter.betas == betas

# Test 3 - Passed
def test_solve_fourier_tf(signal_under_test):
    """Test solve_fourier_tf method with signals from the fixture."""
    uid = 1
    tf = lambda s: s ** 2
    ir = None
    chars = {}
    betas = []
    abstract_filter = AbstractFilter(uid, tf, ir, chars, betas)

    lendata = signal_under_test.length
    fftdata = signal_under_test.mode_f
    freqs = signal_under_test.freqstamps

    print(f"lendata: {lendata}")
    print(f"fftdata: {fftdata}")
    print(f"freqs: {freqs}")

    result = abstract_filter.solve_fourier_tf(lendata, fftdata, freqs)

    print(f"result: {result}")

    assert len(result) == lendata
    assert isinstance(result, np.ndarray)
    # check if the result is not empty.
    assert result.size > 0

# Test 4 - Passed
def test_solve_convolve_ir_exponential():
    """Test solve_convolve_ir with an exponential impulse response."""
    ir = lambda t: np.exp(-t) if t >= 0 else 0
    # Create 100 evenly spaced timestamps from 0 to 5.
    timestamps = np.linspace(0, 5, 100)
    data = np.sin(timestamps)
    filter_obj = AbstractFilter("test", None, ir, {}, [])
    result = filter_obj.solve_convolve_ir(data, timestamps)
    # Assert that the output signal has the same length as the input signal.
    assert len(result) == 100
    # Assert that all elements of the output signal are finite numbers (not NaN or infinity).
    assert np.all(np.isfinite(result))

# Test 5 - Passed
def test_solve_convolve_ir_zero_ir():
    """Test solve_convolve_ir with a zero impulse response."""
    ir = lambda t: 0  # Define an impulse response function that always returns 0.
    timestamps = np.linspace(0, 5, 100)
    # Create a random input signal of 100 values.
    data = np.random.rand(100)
    filter_obj = AbstractFilter("test", None, ir, {}, [])
    result = filter_obj.solve_convolve_ir(data, timestamps)
    # Assert that all elements of the output signal are 0 (since the impulse response is 0).
    assert np.all(result == 0)

# Test 6 - Passed
def test_approx_func_tf_from_ir_exceptions():
    """Test _approx_func_tf_from_ir raises exceptions for invalid input."""
    ir = lambda t: np.exp(-t) if t >= 0 else 0
    filter_obj = AbstractFilter("test", None, ir, {}, [])
    with pytest.raises(Exception, match='no input \:\('):
        filter_obj._approx_func_tf_from_ir(ir, 0, fs=None, time_span=None)
    with pytest.raises(Exception, match='only provide num_samples and one of fs or time_span'):
        filter_obj._approx_func_tf_from_ir(ir, 0, fs=10, time_span=10)
    with pytest.raises(Exception, match='get the scipy error for invalid normalization mode'):
        filter_obj._approx_func_tf_from_ir(ir, 0, fs=10, norm='invalid')

# Test 7 - passed
def test_arbitrary_conflicting_inputs():
    """Test error on dual TF/IR initialization."""
    with pytest.raises(Exception):
        Arbitrary(uid=3, tf=lambda s: 1, ir=lambda t: 1)

# Test 8 - passed
def test_approx_ir_from_tf_convergence():
    """Test _approx_ir_from_tf for convergence with increasing num_samples."""
    filter_instance = AbstractFilter(1, None, None, None, None)
    tf = lambda f: 1 / (1 + 2j * np.pi * f + (2j * np.pi * f)**2)
    t = 1.0
    approx_ir_low = filter_instance._approx_func_ir_from_tf(tf, t, num_samples=500)
    approx_ir_high = filter_instance._approx_func_ir_from_tf(tf, t, num_samples=1000)
    assert np.isclose(approx_ir_low, approx_ir_high, rtol=1e-2) #check that the values are close, indicating convergence.

#Test 9 - failed (complicated error)
def test_approx_ir_from_tf_normalization_modes():
    """Test _approx_ir_from_tf with different normalization modes."""
    filter_instance = AbstractFilter(1, None, None, None, None)
    tf = lambda f: 1.0 + 0j  # Constant transfer function, now complex
    t = 1
    num_samples = 102 # Changed to 102, to avoid num_samples = 1
    backward = filter_instance._approx_func_ir_from_tf(tf, t, norm='backward', num_samples=num_samples)
    ortho = filter_instance._approx_func_ir_from_tf(tf, t, norm='ortho', num_samples=num_samples)
    forward = filter_instance._approx_func_ir_from_tf(tf, t, norm='forward', num_samples=num_samples)
    assert np.isclose(backward, 1.0, rtol=1e-2)
    assert np.isclose(ortho, np.sqrt(num_samples), rtol=1e-2)
    assert np.isclose(forward, num_samples, rtol=1e-2)
    assert backward is not None
    assert ortho is not None
    assert forward is not None

#Test 10 - passed
def test_approx_func_tf_from_ir_approximation():
    """Test _approx_func_tf_from_ir and show it's an approximation."""
    uid = "filter_004"
    tf = None
    ir = lambda t: t + 1
    chars = {}
    betas = []
    af = AbstractFilter(uid, tf, ir, chars, betas)

    f = 1
    fs_low = 5  # Lower sampling frequency
    fs_high = 10  # Higher sampling frequency

    result_low = af._approx_func_tf_from_ir(ir, f, norm="backward", fs=fs_low)
    result_high = af._approx_func_tf_from_ir(ir, f, norm="backward", fs=fs_high)

    print(f"Result with fs={fs_low}: {result_low}")
    print(f"Result with fs={fs_high}: {result_high}")

    # Check if the results are different, indicating approximation
    assert not np.isclose(result_low, result_high, rtol=1e-2)

# Test 11 - Passed
def test_chars2consts_reversible():
    """Tests that chars2consts and sharpfiltercharacteristics (as an approximation) are reversible."""
    characteristics = {'Bpeak': 1.5, 'Nbeta': 11.1, 'phiaccum': 3.5}
    # Call chars2consts and sharpfiltercharacteristics from helpers file
    params = chars2consts(characteristics)
    reconstructed_characteristics = sharpfiltercharacteristics(params['Ap'],params['bp'],params['Bu'])
    # Comparing the values
    assert np.isclose(characteristics['Bpeak'], reconstructed_characteristics['Bpeak'])
    assert np.isclose(characteristics['Nbeta'], reconstructed_characteristics['Nbeta'])
    assert np.isclose(characteristics['phiaccum'], reconstructed_characteristics['phiaccum'])

# Test 12 - Passed
def test_filter_solve_different_methods():
    """Tests Filter solve with different methods."""
    filter = Filter(Ap=0.1, bp=1, Bu=3)
    signal = Signal.from_function(func=np.sin, fs=10, num_samples=100)
    solved_using_tf = filter.solve(signal, method='tf')
    solved_using_ir = filter.solve(signal, method='ir')
    solved_using_ode = filter.solve(signal, method='ode')
    solved_using_fde = filter.solve(signal, method='fde')
    # Comparing the values lengths
    assert len(solved_using_tf) == len(signal)
    assert len(solved_using_ir) == len(signal)
    assert len(solved_using_ode) == len(signal)
    assert len(solved_using_fde) == len(signal)

# Test 13 - 1 passed, 1 failed
def test_arbitrary_init_exceptions():
    """Test exceptions during Arbitrary filter initialization."""
    with pytest.raises(Exception, match="Only one of tf and ir should be used to initialize Filter"):
        Arbitrary(1, tf=lambda s: 1/(1+s), ir=lambda t: np.exp(-t)) # pass (run it alone without the second case)

    with pytest.raises(Exception, match="At least one of tf and ir should be used to initialize Filter"):
        Arbitrary(1)  # Neither tf nor ir provided - should do exception handeling exception (fail)

# Test 14 - passed
def test_rational_initialization_with_coeffs():
    """Test Rational filter initialization with coefficients."""
    uid = 5
    numer = [1, 2, 3]
    denom = [1, 4, 6]
    rational_filter = Rational(uid, coeffs=(numer, denom))
    assert rational_filter.coeffs == (numer, denom) # checking if coefficients are correctly stored
    assert rational_filter.order == 2 # checking if the order is correctly calculated
    assert len(rational_filter.roots[0]) == len(numer) - 1
    assert len(rational_filter.roots[1]) == len(denom) - 1

# Test 15 - passed
def test_rational_initialization_with_roots():
    """Test Rational filter initialization with roots."""
    uid = 2
    zeros = [1, 2]
    poles = [3, 4]
    rational_filter = Rational(uid, roots=(zeros, poles))

    #printing to check the output
    print(f"Rational Filter Roots: {rational_filter.roots}")
    print(f"Rational Filter Order: {rational_filter.order}")
    print(f"Numerator Coefficients Length: {len(rational_filter.coeffs[0])}")
    print(f"Denominator Coefficients Length: {len(rational_filter.coeffs[1])}")

    # assertions comparing
    assert rational_filter.roots == (zeros, poles)
    assert rational_filter.order == 2
    assert len(rational_filter.coeffs[0]) == len(zeros) + 1
    assert len(rational_filter.coeffs[1]) == len(poles) + 1

# Test 16 - passed
def test_rational_invalid_initialization():
    """Test Rational filter initialization with invalid parameters."""
    uid = 4
    with pytest.raises(Exception, match='Invalid rational transfer function filter initialization'):
        Rational(uid) # Should raise an exception as neither coeffs nor roots are given.

# Test 17 - failed
def test_rational_given_coeffs_calculates_roots():
    """Test if Rational filter calculates roots correctly from coefficients."""
    uid = 8
    numer = [1, -3, 2] # Numerator coefficients
    denom = [1, -2, 1] # Denominator coefficients
    rational_filter = Rational(uid, coeffs=(numer, denom))
    calculated_zeros = np.roots(numer).tolist() # Calculate zeros
    calculated_poles = np.roots(denom).tolist() # Calculate poles
    assert np.allclose(rational_filter.roots[0], calculated_zeros) # Check if zeros match calculated zeros.
    assert np.allclose(rational_filter.roots[1], calculated_poles) # Check if poles match calculated poles.

# Test 18 - passed
def test_rational_given_roots_calculates_coeffs():
    """Test if Rational filter calculates coefficients correctly from roots."""
    uid = 2
    zeros = [1, 2] # Zeros of the transfer function
    poles = [3, 4] # Poles of the transfer function
    rational_filter = Rational(uid, roots=(zeros, poles))
    calculated_numer_coeffs = np.polynomial.polynomial.polyfromroots(zeros).tolist() # Calculate numer coeffs from zeros
    calculated_denom_coeffs = np.polynomial.polynomial.polyfromroots(poles).tolist() # Calculate denom coeffs from poles
    assert np.allclose(rational_filter.coeffs[0], calculated_numer_coeffs) # Check numer coeffs match calculated.
    assert np.allclose(rational_filter.coeffs[1], calculated_denom_coeffs) # Check denom coeffs match calculated.

# Test 19 - failed (needs exception handling because it exposes a lot of system details)
def test_rational_solve_diffeq_handles_empty_input():
    """Test if solve_diffeq handles empty input correctly."""
    uid = 7
    numer = [1, 1]
    denom = [1, 2, 1]
    rational_filter = Rational(uid, coeffs=(numer, denom))
    input_data = np.array([]) # Empty input data
    result = rational_filter.solve_diffeq(input_data, fs=1) # Solve difference equation
    assert len(result) == 0 # Check output is empty when input is empty.