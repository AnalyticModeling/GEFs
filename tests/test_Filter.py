import pytest
import numpy as np
from Filter import Filter
from FilterType import Parameterized, Rational, Arbitrary
from Signal import Signal

#======================
# Filters Fixture
#======================

@pytest.fixture(params=[
    # Filters parameters that will be used to create/initialize the filters. (can be extended to include more filters)

    # --- Arbitrary Filters  ---
    # Impulse Response
    {'type': 'impulse_response', 'ir': lambda t: np.exp(-t)},
    # 1.2 Transfer Function
    {'type': 'transfer_function', 'tf': lambda s: 1 / (1 + s + s ** 2)},

    # --- Rational Filters  ---
    # Coefficients
    {'type': 'coefficients', 'coeffs': [[1], [1, 0.5, 0.25]]},
    # Roots
    {'type': 'roots', 'roots': ([-1], [0.5 + 0.5j, 0.5 - 0.5j])},

    # --- Parameterized Filters ---
    # Using Parameters: Ap, bp, Bu, cf.
    {'type': 'parameterized', 'Ap': 0.1, 'bp': 1, 'Bu': 3, 'cf': 10.0},

    # Using Normalized Characteristics
    # 1- phiaccum + Nbeta
    {'type': 'normalized_characteristics', 'phiaccum': np.pi / 2, 'Nbeta': 5.0, 'cf': 1.0},
    # 2- phiaccum + Qerb
    {'type': 'normalized_characteristics', 'phiaccum': np.pi / 2, 'Qerb': 2.0, 'cf': 1.0},
    # 3- phiaccum + Qn
    {'type': 'normalized_characteristics', 'phiaccum': np.pi / 2, 'Qn': 5.0, 'n': 10, 'cf': 1.0},
    # 4- phiaccum + Sbeta
    {'type': 'normalized_characteristics', 'phiaccum': np.pi / 2, 'Sbeta': 15.0, 'cf': 1.0},
    # 5- Nbeta + Qerb
    {'type': 'normalized_characteristics', 'Nbeta': 5.0, 'Qerb': 2.0, 'cf': 1.0},
    # 6- Qn + Qn2 (requires both)
    {'type': 'normalized_characteristics', 'Qn': 5.0, 'Qn2': 0.5, 'cf': 1.0, 'n': 10, 'n2': 20},
    # 7- Sbeta + Qn
    {'type': 'normalized_characteristics', 'Sbeta': 15.0, 'Qn': 5.0, 'n': 10, 'cf': 1.0},

    # Using Unnormalized Characteristics (require cf)
    # 1- Sf + Nf
    {'type': 'unnormalized_characteristics', 'Sf': 10.0, 'Nf': 2.0, 'cf': 2.0},

    # 8- Mixed Characteristics (Normalized and unnormalized)
    { 'type': 'mixed_characteristics','Nbeta': 2.0, 'Sf': 10.0,'cf': 2.0, 'n': 10},

    # --- Multiband Filters ---
    # Multiband with Parameters
    {'type': 'multiband_Parameters', 'Ap': 0.1, 'bp': [0.5, 1, 1.5], 'Bu': [3, 5, 7], 'peak_magndb': 1},

    # Multiband with Characteristics
    {'type': 'multiband_characteristics', 'Bpeak': [0.5, 1, 1.5], 'Nbeta': [15, 10, 5], 'phiaccum': 3.5},

])
def filter_under_test(request):
    """ Create/initialize various Filter objects for testing by requesting the parameters from the filters fixture."""
    params = request.param

    if params['type'] == 'impulse_response':
        return Filter(ir=params['ir'])
    elif params['type'] == 'transfer_function':
        return Filter(tf=params['tf'])
    elif params['type'] == 'coefficients':
        return Filter(coeffs=params['coeffs'])
    elif params['type'] == 'roots':
        return Filter(roots=params['roots'])
    elif params['type'] == 'parameterized':
        return Filter(Ap=params['Ap'], bp=params['bp'], Bu=params['Bu'], cf=params.get('cf', 1.0))
    elif params['type'] == 'normalized_characteristics':
        allowed_chars = {'phiaccum', 'Nbeta', 'Qerb', 'Qn', 'Qn2', 'Sbeta'}
        char_params = {k: v for k, v in params.items() if k in allowed_chars and v is not None}
        n = params.get('n', 10)
        n2 = params.get('n2', 3)
        if len(char_params) != 2:
            pytest.skip(f"Invalid normalized characteristic combination: {char_params.keys()}")
        return Filter(**char_params, cf=params.get('cf', 1.0), n=n, n2=n2)
    elif params['type'] == 'unnormalized_characteristics':
        allowed_chars = {'Sf', 'fpeak', 'Nf', 'ERBf', 'BWndBf', 'BWn2dBf'}
        char_params = {k: v for k, v in params.items() if k in allowed_chars and v is not None}
        if len(char_params) != 2:
            pytest.skip(f"Invalid unnormalized characteristic combination: {char_params.keys()}")
        return Filter(**char_params, cf=params.get('cf', 1.0))
    elif params['type'] == 'mixed_characteristics':
        allowed_chars = {'phiaccum', 'Nbeta', 'Qerb', 'Qn', 'Qn2', 'Sbeta', 'Sf', 'fpeak', 'Nf', 'ERBf', 'BWndBf',
                         'BWn2dBf'}
        char_params = {k: v for k, v in params.items() if k in allowed_chars and v is not None}
        n = params.get('n', 10)
        n2 = params.get('n2', 3)
        if len(char_params) != 2:
            pytest.skip(f"Invalid mixed characteristic combination: {char_params.keys()}")
        return Filter(**char_params, cf=params.get('cf', 1.0), n=n, n2=n2)
    elif params['type'] == 'multiband_Parameters':
        return Filter.multiband_consts(Ap=params['Ap'], bp=params['bp'], Bu=params['Bu'],
                                           peak_magndb=params['peak_magndb'])
    elif params['type'] == 'multiband_characteristics':
        return Filter.multiband_chars(Bpeak=params['Bpeak'], Nbeta=params['Nbeta'], phiaccum=params['phiaccum'])
    return None

#======================
# Signals Fixture
#======================

@pytest.fixture(params=[
    # Signals parameters that will be used to create/initialize the signals. (can be extended to include more signals)

    # Time domain signals
    {'type': 'time', 'mode': 't', 'length': 100, 'fs': 1.0},
    # Description: Standard time-domain signal with samples at regular intervals.

    {'type': 't_tilde', 'mode': 'ttilde', 'length': 500, 'fs': 10.0},
    # Description: Time-domain signal with a modified time axis or representation.

    # Frequency domain signals
    {'type': 'frequency', 'mode': 'f', 'length': 64, 'fs': 2.0, 'evenlen': True},
    # Description: Signal in the frequency domain (DFT) with linear frequency axis.

    {'type': 'angular', 'mode': 'w', 'length': 128, 'fs': 4.0, 'evenlen': False},
    # Description: Signal in the frequency domain with angular frequency (rad/s) axis.

    # Factory method generated signals
    {'type': 'function_generated', 'factory': 'from_function',
     'func': lambda t: np.sin(2 * np.pi * 5 * t), 'fs': 10.0, 'num_samples': 100},
    # Description: Signal generated by evaluating a user-defined function.

    {'type': 'linear_chirp', 'factory': 'linear_chirp','f_init': 1.0, 'f_final': 10.0, 'fs': 100.0, 'num_samples': 200},
    # Description: Signal with linearly increasing/decreasing frequency over time.

    {'type': 'instantaneous_frequency_generated', 'factory': 'from_instantaneous_frequency',
     'freq_func': lambda t: t / 10, 'fs': 10.0, 'num_samples': 50},
    # Description: Signal generated by specifying its instantaneous frequency as a function of time.

    # Boundary cases
    {'type': 'minimal_length', 'mode': 't', 'length': 1, 'fs': 1.0},
    # Description: Signal with the minimum possible length (1 sample).

    {'type': 'nyquist_frequency_limit', 'mode': 't', 'data': [np.sin(2 * np.pi * 0.45 * 1.0 * t) for t in range(100)], 'fs': 1.0},
    # Description: Signal with a frequency close to the Nyquist limit.

    {'type': 'zero_valued', 'mode': 'f', 'data': np.zeros(10), 'fs': 1.0},
    # Description: Signal where all data points are zero.

])
def signal_under_test(request):
    """ Create/initialize various Signals objects for testing by requesting the parameters from the signals fixture."""
    param = request.param

    if 'data' in param:
        return Signal(mode=param['mode'], data=param['data'], fs=param['fs'])

    elif 'factory' in param:
        factory_method = getattr(Signal, param['factory'])
        args = {k: v for k, v in param.items() if k not in ['type', 'factory']}
        return factory_method(**args)

    # Handles signals that have a length, and need random data generated.
    elif 'length' in param:
        return Signal(
            mode=param['mode'],
            data=np.random.randn(param['length']),
            fs=param['fs'],
            evenlen=param.get('evenlen', True))

    else:
        raise ValueError(f"Invalid signal parameter set: {param}")

#======================
# Core Tests
#======================

# Test 1 - Passed
def test_filter_general_properties(filter_under_test):
    """Validates fundamental properties of the Filter objects."""
    # Check for some general attributes (can be extended)
    assert hasattr(filter_under_test, 'uid'), "Filter object should have a 'uid' attribute."
    assert isinstance(filter_under_test.uid, int), "UID should be an integer."
    assert hasattr(filter_under_test, 'cf'), "Filter object should have a 'cf' attribute."
    assert isinstance(filter_under_test.cf, (int, float)), "cf should be a number."
    assert filter_under_test.cf >= 0, "cf should be non-negative."
    if hasattr(filter_under_test, 'cf'):
        assert filter_under_test.cf >= 0, "cf should be non-negative."
    if hasattr(filter_under_test, 'Bpeak'):
        assert filter_under_test.Bpeak > 0, "Bpeak should be positive."

# Test 2 - Passed
def test_invalid_filter_initialization():
    """Test various invalid initialization scenarios"""
    # Conflicting specifications
    with pytest.raises(Exception):
        Filter(ir=lambda t: t, coeffs=[[1], [1]])

    # Missing required parameters
    with pytest.raises(Exception):
        Filter(fpeak=1.0)

    # Invalid filter type
    with pytest.raises(Exception):
        Filter(type='invalid', Ap=0.1, bp=1.0, Bu=3.0)

    # empty filter type
    with pytest.raises(Exception):
        Filter(type='', Ap=0.1, bp=1.0, Bu=50)

    # high values filter attributes
    filter=Filter( Ap=100000000000, bp=500.767, Bu=-10000000)
    assert isinstance(filter, Filter)

    with pytest.raises(Exception):
        # Conflicting specifications
        Filter(tf=lambda s: 1, ir=lambda t: 1)

    with pytest.raises(Exception):
        # Incomplete parameterized filter
        Filter(Ap=0.1, bp=1.0)  # Missing Bu

    with pytest.raises(Exception):
        # Invalid characteristic frequency
        Filter(Nf=10, Qerb=5)  # Missing cf

# Test 3 - Passed (But it include the mixed Characteristics (Normalized and unnormalized))
def test_get_computed_normalized_chars(filter_under_test):
    """get computed characteristics based on provided normalized characteristics."""
    if (filter_under_test is not None and filter_under_test.in_terms_of_normalized is True and isinstance(filter_under_test.filter, Parameterized)):
        characteristics = filter_under_test.get_computed_chars()
        assert isinstance(characteristics, dict)
        # Print computed characteristics
        print("\n\nComputed Normalized Characteristics:")
        for char_name, char_value in characteristics.items():
            print(f"{char_name}: {char_value}")

# Test 4 - Passed (But it include the mixed Characteristics (Normalized and unnormalized))
def test_get_computed_unnormalized_chars(filter_under_test):
    """get computed characteristics based on provided unnormalized characteristics."""
    if filter_under_test is not None and filter_under_test.in_terms_of_normalized is False and isinstance(filter_under_test.filter, Parameterized):
        chars = filter_under_test.get_computed_unnormalized_chars()
        assert isinstance(chars, dict)
        # Print computed characteristics
        print("\n\nComputed Unnormalized Characteristics:")
        for char_name, char_value in chars.items():
            print(f"{char_name}: {char_value}")

# Test 5 - Passed ( but it also get some of computed characteristics from the filter that initialized given parameters )
def test_filter_get_orig_chars(filter_under_test):
    """Test get_orig_chars for Parameterized filters, ensuring it only returns the original characteristics."""
    if filter_under_test is not None and isinstance(filter_under_test.filter, Parameterized):
        # Get the original characteristics from the filter's orig_chars attribute
        orig_chars = filter_under_test.filter.orig_chars

        # Get the characteristics returned by get_orig_chars()
        returned_chars = filter_under_test.get_orig_chars()

        # Assert that the returned characteristics are the same as the original ones (keys)
        assert set(returned_chars.keys()) == set(orig_chars.keys())

        # Assert that the values are the same
        for key in orig_chars:
            assert returned_chars[key] == orig_chars[key]

        # Print the original characteristics
        print("\n\nOriginal Characteristics:")
        for char_name, char_value in returned_chars.items():
            print(f"{char_name}: {char_value}")

# Test 6 - Passed ( but it also get some of computed parameters from the filter that initialized given characteristics )
def test_filter_get_consts(filter_under_test):
    """Test get_consts for Parameterized filters, ensuring it only returns the original constants."""
    if filter_under_test is not None and isinstance(filter_under_test.filter, Parameterized) :
        # Get the original parameters from the filter's consts attribute
        orig_consts = filter_under_test.filter.params

        # Get the parameters returned by get_consts()
        returned_consts = filter_under_test.get_consts()

        # Assert that the returned parameters are the same as the original ones (keys)
        assert set(returned_consts.keys()) == set(orig_consts.keys())

        # Assert that the values are the same
        for key in orig_consts:
            assert returned_consts[key] == orig_consts[key]

        # Print the original parameters
        print("\n\nOriginal Parameters:")
        for param_name, param_value in returned_consts.items():
            print(f"{param_name}: {param_value}")

# Test 7 - 108 passed, 20 failed and remaining facing performance problem.
def test_filters_solveing_signals(filter_under_test, signal_under_test):
    """Test the solve method with various all possible combinations between filters and signals using all methods."""
    if filter_under_test is not None and signal_under_test is not None:
        # Test with default method (None)
        output = filter_under_test.solve(signal_under_test)
        assert isinstance(output, Signal)
        assert len(output) == len(signal_under_test)
        assert output.fs == pytest.approx(signal_under_test.fs, rel=1e-9)

        # Test with 'tf' method
        output = filter_under_test.solve(signal_under_test, method='tf')
        assert isinstance(output, Signal)
        assert len(output) == len(signal_under_test)
        assert output.fs == pytest.approx(signal_under_test.fs, rel=1e-9)

        # Test with 'ir' method if available
        if hasattr(filter_under_test.filter, 'ir'):
            output = filter_under_test.solve(signal_under_test, method='ir')
            assert isinstance(output, Signal)
            assert len(output) == len(signal_under_test)
            assert output.fs == pytest.approx(signal_under_test.fs, rel=1e-9)

        # Test with 'ode' method if available
        if isinstance(filter_under_test.filter, (Rational, Parameterized)):
            try:
                output = filter_under_test.solve(signal_under_test, method='ode')
                assert isinstance(output, Signal)
                assert len(output) == len(signal_under_test)
                assert output.fs == pytest.approx(signal_under_test.fs, rel=1e-9)
            except Exception as e:
                # If 'ode' method is not applicable, catch the exception and print a warning
                print(f"Warning: 'ode' method not applicable - {e}")

        # Test with 'fde' method if available
        if isinstance(filter_under_test.filter, Parameterized):
            try:
                output = filter_under_test.solve(signal_under_test, method='fde')
                assert isinstance(output, Signal)
                assert len(output) == len(signal_under_test)
                assert output.fs == pytest.approx(signal_under_test.fs, rel=1e-9)
            except Exception as e:
                # If 'fde' method is not applicable, catch the exception and print a warning
                print(f"Warning: 'fde' method not applicable - {e}")

        # Test with non-Signal input (provide fs)
        if not isinstance(signal_under_test, Signal):
            output = filter_under_test.solve(signal_under_test.mode_t, fs=signal_under_test.fs)
            assert isinstance(output, Signal)
            assert len(output) == len(signal_under_test)
            assert output.fs == pytest.approx(signal_under_test.fs, rel=1e-9)

        # Test with invalid method
        with pytest.raises(Exception) as exception:
            filter_under_test.solve(signal_under_test, method='invalid')
        assert "Invalid solution method" in str(exception.value)

# Test 8 - Passed
def test_bode_plot_properties(filter_under_test):
    """Test bode_plot method properties"""
    if filter_under_test is not None:
        # Call the plot with show=False to prevent the plot from being displayed.
        x_axis, magnitudes, phases = filter_under_test.bode_plot(show=False)
        # Assert that the returned values are of the correct types and lengths
        assert isinstance(x_axis, np.ndarray)
        assert isinstance(magnitudes, np.ndarray)
        assert isinstance(phases, np.ndarray)
        assert len(x_axis) == len(magnitudes) == len(phases)

# Test 9 - Passed
def test_frequency_real_imaginary_plot_properties(filter_under_test):
    """Test frequency_real_imag_plot properties."""
    if filter_under_test is not None:
        """
        The method is expected to return:
          x_axis: The x-axis values of the plot, which are typically frequencies.
          reals: The real parts of the frequency response.
          imaginary: The imaginary parts of the frequency response.
        """
        # Call the plot with show=False to prevent the plot from being displayed.
        x_axis, reals, imaginary = filter_under_test.frequency_real_imag_plot(show=False)
        # Assert that the returned values are of the correct types and lengths
        assert isinstance(x_axis, np.ndarray)
        assert isinstance(reals, list)
        assert isinstance(imaginary, list)
        assert len(x_axis) == len(reals) == len(imaginary)

# Test 10 - Passed
def test_nichols_plot_properties(filter_under_test):
    """Test nichols_plot properties."""
    if filter_under_test is not None:
        """
        The method is expected to return:
          x_axis: The x-axis values of the Nichols plot, which are typically frequencies.
          phases: The phases (in degrees or radians) of the frequency response.
          magnitudes: The magnitudes (in dB) of the frequency response.
        """
        # Call the plot with show=False to prevent the plot from being displayed.
        x_axis, phases, magnitudes = filter_under_test.nichols_plot(show=False)
        # Assert that the returned values are of the correct types and lengths
        assert isinstance(x_axis, np.ndarray)
        assert isinstance(phases, np.ndarray)
        assert isinstance(magnitudes, np.ndarray)
        assert len(x_axis) == len(phases) == len(magnitudes)

# Test 11 - Passed
def test_nyquist_plot_properties(filter_under_test):
    """Test nyquist_plot properties."""
    if filter_under_test is not None:
        """
        The method is expected to return:
          x_axis: The x-axis values of the Nyquist plot, which are typically frequencies.
          reals: The real parts of the frequency response.
          imaginary: The imaginary parts of the frequency response.
        """
        # Call the plot with show=False to prevent the plot from being displayed.
        x_axis, reals, imaginary = filter_under_test.nyquist_plot(show=False)
        # Assert that the returned values are of the correct types and lengths
        assert isinstance(x_axis, np.ndarray)
        assert isinstance(reals, np.ndarray)
        assert isinstance(imaginary, np.ndarray)
        assert len(x_axis) == len(reals) == len(imaginary)

# Test 12 - 14 passed, 2 failed and cause an overflow without exception
def test_impulse_response_plot_properties(filter_under_test):
    """Test impulse_response_plot properties."""
    if filter_under_test is not None:
        """
        The method is expected to return:
          x_axis: The x-axis values of the plot, which are typically time values.
          response: The impulse response values.
        """
        # Call the plot with show=False to prevent the plot from being displayed.
        x_axis, response = filter_under_test.impulse_response_plot(show=False)
        # Assert that the returned values are of the correct types and lengths
        assert isinstance(x_axis, np.ndarray)
        assert isinstance(response, np.ndarray)
        assert len(x_axis) == len(response)

# Test 13 - Passed
def test_pole_zero_plot_properties(filter_under_test):
    """Test pole_zero_plot properties."""
    if filter_under_test is not None and isinstance(filter_under_test.filter, (Rational, Parameterized)):
        """
        The method is expected to return:
          zeros: The zeros of the transfer function.
          poles: The poles of the transfer function.
        """
        # Call the plot with show=False to prevent the plot from being displayed.
        zeros, poles = filter_under_test.pole_zero_plot(show=False)
        # Assert that the returned values are of the correct types
        assert isinstance(zeros, list)
        assert isinstance(poles, list)
    elif filter_under_test is not None and not isinstance(filter_under_test.filter, (Rational, Parameterized)):
        # For other filter types, assert that an exception is raised
        with pytest.raises(Exception):
            filter_under_test.pole_zero_plot(show=False)

# Test 14 - Passed
def test_Qn_plot_properties(filter_under_test):
    """Test Qn_plot properties."""
    if filter_under_test is not None and isinstance(filter_under_test.filter, Parameterized):
        """
        Call Qn_plot with show=False to prevent the plot from being displayed.
        The method is expected to return:
          x_axis: The x-axis values of the plot, which are typically the 'n' values (in dB).
          Qns: The corresponding Qn values.
        """
        # Call the plot with show=False to prevent the plot from being displayed.
        x_axis, Qns = filter_under_test.Qn_plot(show=False)
        # Assert that the returned values are of the correct types and lengths
        assert isinstance(Qns, list)
        assert len(x_axis) == len(Qns)

    elif filter_under_test is not None and not isinstance(filter_under_test.filter, Parameterized):
        # For other filter types, assert that an exception is raised
        with pytest.raises(Exception, match='Qn plot not guaranteed' ):
            filter_under_test.Qn_plot(show=False)

# Test 15 - failed (could not solve or handel solving and empty signal data)
def test_edge_cases_empty_signal_data():
    """Test boundary conditions and edge cass"""
    # Empty signal
    empty_signal = Signal(mode='t', data=[], fs=1.0)
    filter = Filter(coeffs=[[1], [1]])
    filter.solve(empty_signal)

# Test 16 - passed
def test_parameterized_filter_characteristics():
    """Test parameter conversion from characteristics to parameters"""
    # Known characteristic -> parameter relationship
    f = Filter(Nbeta=10, Qerb=20, cf=1.0)
    consts = f.get_consts()

    # Validate parameter relationships
    assert np.isclose(consts['bp'], 1.0)
    assert consts['Ap'] > 0
    assert consts['Bu'] > 0

    # Verify computed characteristics match requested
    chars = f.get_computed_chars()
    assert np.isclose(chars['Nbeta'], 10, rtol=0.1)
    assert np.isclose(chars['Qerb'], 20, rtol=0.1)

# Test 17 - passed
def test_characteristic_error_validation_for_non_parameterized():
    """Test error conditions and edge cases for characteristic_error()."""
    # Non-Parameterized filters should raise exception
    with pytest.raises(Exception):
        Filter(coeffs=[[1], [1, 1]]).characteristic_error()


# Test 18 - failed (parameterized without original chars (initialized with consts) should raise exception)
def test_characteristic_error_validation_for_parameterized_using_parms():
    """Test error conditions and edge cases for characteristic_error()."""
    with pytest.raises(Exception):
        Filter(Ap=0.1, bp=1, Bu=3).characteristic_error()

# Test 19 - passed
def test_characteristic_error_calculation():
    """Validate error calculation accuracy with known values."""
    # Create filter with known characteristics
    filter = Filter(Nbeta=5.0, Qerb=20.0, cf=1.0)

    # Get error data
    error_data = filter.characteristic_error(show=False)
    original_characteristics = filter.get_orig_chars()
    computed_characteristics = filter.get_computed_chars()

    # Verify error calculation
    for characteristic in original_characteristics:
        expected = (computed_characteristics[characteristic] / original_characteristics[characteristic]) - 1
        #if the difference between the two numbers is 0.01 or less (the number can be modified).
        assert np.isclose(error_data[characteristic], expected, atol=0.01)

# Test 20 - passed
def test_characteristic_error_attribute_appearance():
    """Test which attribute are included to calculate the characteristic error"""
    filter = Filter(Sf=10.0, Nf=2.0, cf=2.0)

    error_data = filter.characteristic_error(show=False)

    # Original chars would be Sbeta and Nbeta internally
    assert 'Sf' in error_data
    assert 'Nf' in error_data
    assert 'Sbeta' not in error_data
    assert 'Nbeta' not in error_data

# Test 21 - failed (it should raise exception because the characteristic_error feature not implemented for multiband characteristic filters)
def test_multiband_characteristic_error():
    """Test error handling for multiband characteristic filters."""
    filt = Filter.multiband_chars(
        Bpeak=[0.5, 1, 1.5],
        Nbeta=[15, 10, 5],
        phiaccum=3.5
    )
    with pytest.raises(NotImplementedError):
        filt.characteristic_error()