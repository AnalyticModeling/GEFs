import pytest
import numpy as np
import scipy.fft
from src.GEFs_core.Signal import Signal

#======================
# Signal Fixture
#======================
# (can be extended to include more signals)
@pytest.fixture
def sample_signals():
    """Fixture providing sample Signal objects for testing."""
    # sampling frequency
    fs = 1000
    t = np.linspace(0, 1, fs, endpoint=False)
    data = np.sin(2 * np.pi * 5 * t)
    return Signal(mode='t', data=data, fs=fs)

#======================
# Code Tests
#======================

# Test 1 - passed
def test_signal_initialization_defaults():
    """Test default initialization of Signal."""
    signal = Signal()
    assert signal.fs == 1
    assert len(signal) == 9
    assert np.allclose(signal.get_data('t'), [0] * 9)

# Test 2 - passed
def test_signal_initialization_time_domain():
    """Test initialization with time-domain data."""
    fs = 10
    data = [1, 2, 3, 4, 5]
    signal = Signal(mode='t', data=data, fs=fs)
    assert signal.fs == fs
    assert len(signal) == len(data)
    assert np.array_equal(signal.get_data('t'), data)

# Test 3 - passed
def test_signal_initialization_frequency_domain_even():
    """Test initialization with frequency-domain data (even length)."""
    data_f = [1, 2, 3]
    fs = 10
    signal = Signal(mode='f', data=data_f, fs=fs, evenlen=True)
    expected_length = 2 * len(data_f) - 2
    assert signal.fs == fs
    assert len(signal) == expected_length
    assert np.array_equal(signal.get_data('f'), data_f)

# Test 4 - passed
def test_signal_from_function_time():
    """Test creating a Signal from a custom function in time-domain."""
    func = lambda t: t ** 2
    fs = 10
    num_samples = 5
    signal = Signal.from_function(mode='t', func=func, fs=fs, num_samples=num_samples)
    t_values = np.arange(num_samples) / fs
    assert np.allclose(signal.get_data('t'), func(t_values))

# Test 5 - passed
def test_signal_from_function_frequency():
    """Test creating a Signal from a custom function in frequency-domain."""
    fs = 10
    num_samples = 5
    func = lambda f: 1 / (1 + f ** 2)
    signal = Signal.from_function(mode='f', func=func, fs=fs, num_samples=num_samples, evenlen=True)
    freq_values = scipy.fft.rfftfreq(2 * num_samples - 2, 1 / fs)
    assert np.allclose(signal.get_data('f'), func(freq_values))

# Test 6 - passed
def test_signal_linear_chirp():
    """Test creating a linear chirp Signal."""
    signal = Signal.linear_chirp(f_init=10, f_final=20, fs=1000, num_samples=100)
    assert signal.fs == 1000
    assert len(signal) == 100
    assert signal.mode == 't'

# Test 7 - passed
def test_linear_chirp_angular_frequency():
    """Test Signal.linear_chirp with angular frequencies."""
    signal = Signal.linear_chirp(w_init=2 * np.pi, w_final=6 * np.pi, fs=10, num_samples=5)
    timestamps = np.arange(5) / 10
    expected_data = [np.cos(np.pi * t * (1 * (2 - t / 0.4) + 3 * (t / 0.4))) for t in timestamps]
    assert np.allclose(signal.get_data('t'), expected_data)

# Test 8 - passed
def test_signal_from_instantaneous_frequency():
    """Test creating a Signal using instantaneous frequency."""
    fs = 10
    freq_func = lambda t: 1 + t
    signal = Signal.from_instantaneous_frequency(freq_func=freq_func, fs=fs, num_samples=5)
    assert signal.fs == fs
    assert len(signal) == 5

# Test 9 - passed
def test_signal_initialization_invalid_mode():
    """Test initialization with an invalid mode."""
    with pytest.raises(Exception):
        Signal(mode='invalid', data=[1, 2, 3])

# Test 10 - failed
def test_signal_resample(sample_signals):
    """Test resampling the Signal."""
    signal = sample_signals
    new_fs = 2000
    resampled_signal = signal.resample(new_fs)
    assert resampled_signal.fs == new_fs #passed
    # resampled one should get more lengths because it had higher fs (the error happened because resampled_signal was 501, and the length of the original signal was 1000)
    assert len(resampled_signal) > len(signal) #failed

# Test 11 - failed ( to make it passed i need to remove 'window=' argument from the original method (spectral_entropy))
def test_signal_spectral_entropy(sample_signals):
    """Test computation of spectral entropy."""
    signal = sample_signals
    entropy = signal.spectral_entropy()
    assert isinstance(entropy, float)
    # Spectral entropy is typically normalized to a range between 0 and 1.
    assert 0 <= entropy <= 1

# Test 12 - passed
def test_signals_addition():
    """Test addition of two Signals."""
    signal1 = Signal(mode='t', data=[1, 2, 3], fs=10)
    signal2 = Signal(mode='t', data=[4, 5, 6], fs=10)
    result = signal1.__add__( signal2)
    expected_data = [5, 7, 9]
    assert np.allclose(result.get_data('t'), expected_data)

# Test 13 - passed
def test_signal_addition_different_lengths():
    """Test Signal addition with different lengths."""
    signal1 = Signal(mode='t', data=[1, 2, 3], fs=10)
    signal2 = Signal(mode='t', data=[4, 5], fs=10)
    result = signal1.__add__( signal2)
    expected_data = [5, 7, 3] #padding with zeros
    assert np.allclose(result.get_data('t'), expected_data)

# Test 14 - passed
def test_signal_addition_different_fs():
    """Test Signal addition with different sampling rates."""
    signal1 = Signal(mode='t', data=[1, 2, 3], fs=10)
    signal2 = Signal(mode='t', data=[4, 5, 6], fs=20)
    with pytest.raises(Exception) as exception_info:
        signal1.__add__( signal2)
    assert str(exception_info.value) == 'Cannot add two signals with different sampling rates'

# Test 15 - passed
def test_signal_scalar_addition():
    """Test Signal right addition."""
    signal = Signal(mode='t', data=[1, 2, 3], fs=10)
    result = signal.__add__( 5)
    expected_data = [6, 7, 8]
    assert np.allclose(result.get_data('t'), expected_data)

# Test 16 - passed
def test_signal_subtraction():
    """Test Signal subtraction."""
    signal1 = Signal(mode='t', data=[1, 2, 3], fs=10)
    signal2 = Signal(mode='t', data=[4, 5, 6], fs=10)
    result = signal1.__sub__( signal2)
    expected_data = [-3, -3, -3]
    assert np.allclose(result.get_data('t'), expected_data)

# Test 17 - passed
def test_signal_multiply_scalar(sample_signals):
    """Test multiplication of a scalar with the Signal."""
    signal = sample_signals
    result = signal.__mul__(2)
    assert np.allclose(result.get_data('t'), signal.get_data('t') * 2)

# Test 18 - passed
def test_signal_multiplication_with_similar_properties_signals():
    """Test Signal multiplication same sampling rates and same data lengths."""
    signal1 = Signal(mode='t', data=[1, 2, 3], fs=10)
    signal2 = Signal(mode='t', data=[4, 5, 6], fs=10)
    result = signal1 * signal2
    expected_data = [4, 10, 18]
    assert np.allclose(result.get_data('t'), expected_data)

# Test 19 - passed
def test_signal_multiply_signals_different_sampling_rates(sample_signals):
    """Test multiplication of two Signals with different sampling rates."""
    signal1 = Signal(mode='t', data=[1, 2, 3, 4, 5], fs=20)
    signal2 = Signal(mode='t', data=[2, 3, 4, 5, 6, 7], fs=10)

    with pytest.raises(Exception, match= "Cannot multiply two signals with different sampling rates"):
        signal1.__mul__(signal2)

# Test 20 - passed
def test_signal_multiply_signals_different_data_lengths():
    """Test multiplication of two Signals that have different data lengths"""
    signal1 = Signal(mode='t', data=[1, 2, 3, 4, 5], fs=10)
    signal2 = Signal(mode='t', data=[2, 3, 4, 5, 6, 7], fs=10)

    result = signal1.__mul__(signal2)
    expected_result = np.array([1 * 2, 2 * 3, 3 * 4, 4 * 5, 5 * 6])  # Truncated to shorter length
    assert np.allclose(result.get_data('t'), expected_result)

# Test 21 - passed
def test_signal_division():
    """Test Signal division."""
    signal1 = Signal(mode='t', data=[4, 10, 18], fs=10)
    signal2 = Signal(mode='t', data=[2, 5, 6], fs=10)
    result = signal1.__truediv__ (signal2)
    expected_data = [2, 2, 3]
    assert np.allclose(result.get_data('t'), expected_data)

# Test 22 - passed
def test_signal_modulo():
    """Test Signal modulo."""
    signal = Signal(mode='t', data=[5, 7, 9], fs=10)
    result = signal.__mod__(3)
    expected_data = [2, 1, 0]
    assert np.allclose(result.get_data('t'), expected_data)

# Test 23 - passed
def test_signal_floor_division():
    """Test Signal floor division."""
    signal = Signal(mode='t', data=[5, 7, 9], fs=10)
    result = signal.__floordiv__(3)
    expected_data = [1, 2, 3]
    assert np.allclose(result.get_data('t'), expected_data)

# Test 24 - passed
def test_signal_power():
    """Test Signal power."""
    signal = Signal(mode='t', data=[1, 2, 3], fs=10)
    result = signal.__pow__(2)
    expected_data = [1, 4, 9]
    assert np.allclose(result.get_data('t'), expected_data)

# Test 25 - passed
def test_signal_negation():
    """Test Signal negation."""
    signal = Signal(mode='t', data=[1, -2, 3], fs=10)
    result = signal.__neg__()
    expected_data = [-1, 2, -3]
    assert np.allclose(result.get_data('t'), expected_data)

# Test 26 - passed
def test_signal_absolute():
    """Test Signal absolute value."""
    signal = Signal(mode='t', data=[1, -2, 3], fs=10)
    result = signal.__abs__()
    expected_data = [1, 2, 3]
    assert np.allclose(result.get_data('t'), expected_data)

# Test 27 - failed ( needs an exception handling )
def test_signal_invalid_resample(sample_signals):
    """Test resampling to an invalid sample rate."""
    signal = sample_signals
    signal.resample(-1)

# Test 28 - passed
def test_pad_time_series_same_length():
    """Test padding when signal1 is longer than signal2."""
    signal1 = Signal(mode='t', data=[1, 2, 3, 4, 5], fs=10)
    signal2 = Signal(mode='t', data=[6, 7, 8], fs=10)

    padded_signal1, padded_signal2 = signal1._pad_time_series_to_same_length(signal2.mode_t)

    assert len(padded_signal1) == len(padded_signal2)  # Check same length
    assert list(padded_signal1) == [1, 2, 3, 4, 5]  # Check signal1 unchanged
    assert list(padded_signal2) == [6, 7, 8, 0, 0]  # Check signal2 padded correctly

# Test 29 - passed
def test_pad_time_series_equal_length():
    """Test when signal1 and signal2 have the same length."""
    signal1 = Signal(mode='t', data=[1, 2, 3], fs=10)
    signal2 = Signal(mode='t', data=[4, 5, 6], fs=10)

    padded_signal1, padded_signal2 = signal1._pad_time_series_to_same_length(signal2.mode_t)

    assert len(padded_signal1) == len(padded_signal2)  # Check same length
    assert list(padded_signal1) == [1, 2, 3]  # Check signal1 unchanged
    assert list(padded_signal2) == [4, 5, 6]  # Check signal2 unchanged

# Test 30 - passed
def test_autocorrelate():
    """Test autocorrelation of Signal."""
    signal = Signal(mode='t', data=[1, 2, 3, 2, 1], fs=10)
    autocorr = signal.autocorrelate()
    expected_autocorr = [19, 16, 10, 4, 1]
    assert np.allclose(autocorr, expected_autocorr)

# Test 31 - passed
def test_autocorrelate_zero_signal():
    """Test autocorrelation of a zero-valued Signal."""
    signal = Signal(mode='t', data=[0, 0, 0], fs=10)
    autocorr = signal.autocorrelate()
    expected_autocorr = [0.0, 0.0, 0.0]
    assert np.allclose(autocorr, expected_autocorr)

# Test 32 - passed
def test_crosscorrelate():
    """Test crosscorrelation of two Signals."""
    signal1 = Signal(mode='t', data=[1, 2, 3], fs=10)
    signal2 = Signal(mode='t', data=[4, 5, 6], fs=10)
    crosscorr = signal1.crosscorrelate(signal2)
    expected_crosscorr = [6, 17, 32, 23, 12]
    assert np.allclose(crosscorr, expected_crosscorr)

# Test 33 - passed
def test_crosscorrelate_different_lengths():
    """Test crosscorrelation of two Signals with different lengths."""
    signal1 = Signal(mode='t', data=[1, 2, 3], fs=10)
    signal2 = Signal(mode='t', data=[4, 5], fs=10)
    crosscorr = signal1.crosscorrelate(signal2)
    expected_crosscorr = [5, 14, 23, 12]
    assert np.allclose(crosscorr, expected_crosscorr)

# Test 34 - passed
def test_crosscorrelate_invalid_input():
    """Test crosscorrelation with invalid input."""
    signal = Signal(mode='t', data=[1, 2, 3], fs=10)
    with pytest.raises(Exception) as exception_info:
        signal.crosscorrelate([1, 2, 3])
    assert str(exception_info.value) == 'can only cross correlate two signals'

# Test 35 - passed
def test_matmul_operator():
    """Test the @ (matmul) operator for crosscorrelation."""
    signal1 = Signal(mode='t', data=[1, 2, 3], fs=10)
    signal2 = Signal(mode='t', data=[4, 5, 6], fs=10)
    crosscorr = signal1.__matmul__(signal2)
    expected_crosscorr = [6, 17, 32, 23, 12]
    assert np.allclose(crosscorr, expected_crosscorr)

# Test 36 - passed
def test_as_sound(tmp_path):
    """Test saving Signal as sound file."""
    signal = Signal(mode='t', data=[100, 200, 300], fs=10)
    filename = tmp_path / "test_sound.wav"
    signal.as_sound(filename)
    assert filename.exists()
    rate, data = scipy.io.wavfile.read(filename)
    assert rate == 10
    assert np.allclose(data, [100, 200, 300])