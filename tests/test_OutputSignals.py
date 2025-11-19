import pytest
import numpy as np
from unittest.mock import patch  # For mocking matplotlib
from src.Signal import Signal
from src.OutputSignals import OutputSignals

#======================
# Signals Fixture
#======================
# (can be extended to include more signals)
@pytest.fixture
def sample_signals():
    """Fixture to create sample Signal objects with data for testing."""
    fs = 1000
    t = np.linspace(0, 1, fs, endpoint=False)
    data_in = np.sin(2 * np.pi * 5 * t)
    data_1 = np.cos(2 * np.pi * 10 * t)
    data_2 = np.sin(2 * np.pi * 15 * t)
    data_3 = np.cos(2 * np.pi * 20 * t)

    # signal_in is the input signal, the root and the original signal before any processing
    signal_in = Signal(data=data_in, fs=fs)
    signal_1 = Signal(data=data_1, fs=fs)
    signal_2 = Signal(data=data_2, fs=fs)
    signal_3 = Signal(data=data_3, fs=fs)

    return signal_in, signal_1, signal_2, signal_3

#======================
# Graph Fixture
#======================
# (can be modified as needed but in way that matches the signals number)
@pytest.fixture
def sample_graph():
    """Fixture to create a sample graph (dictionary) for testing."""
    return type('obj', (object,), {"parent": {1: 0, 2: 0, 3: 2}}) # in the sample graph: 1 and 2 are children of 0, 0 is the root (-1) and 3 is child of 2

#======================
# Core Tests
#======================

# Test 1 - failed (need review)
def test_output_signals_result_and_ids_conversion():
    """Test OutputSignals output as a dictionary and the identification conversion between uid and graphid."""
    signal1 = Signal(mode='t', data=[1, 2, 3], fs=10)
    signal2 = Signal(mode='t', data=[4, 5, 6], fs=10)
    graph_dict = {1: 0, 2: 1}
    output_signals = OutputSignals([signal1, signal2], graph_dict)
    assert output_signals.insignal == signal1
    assert output_signals.outsignals == [signal2]
    assert output_signals.graph == graph_dict
    assert len(output_signals) == 1
    # the 2 follwing assersion causing error, i think it needs a review if i did the logic of _uid2graphid and _graphid2uid correctly or not to expect the output
    assert output_signals._uid2graphid == {-1: 0, 2: 1}
    assert output_signals._graphid2uid == {0: -1, 1: 2}

# Test 2 - passed
def test_outputsignals_init_type_error(sample_signals, sample_graph):
    """Test OutputSignals initialization with non-Signal objects."""
    signal_in, signal_1, signal_2, signal_3 = sample_signals
    with pytest.raises(TypeError) as exception_info:
        OutputSignals([signal_in, signal_1, "not a signal"], sample_graph)
    assert str(exception_info.value) == "All inputs to OutputSignals should be Signal objects"

# Test 3 - failed ( did not raise an exception when duplicated uid was there)
def test_output_signals_initialization_uid_error():
    """Test OutputSignals initialization with non-unique uids."""
    signal1 = Signal(mode='t', data=[1, 2, 3], fs=10)
    signal2 = Signal(mode='t', data=[4, 5, 6], fs=10)
    signal2.uid = 1
    with pytest.raises(Exception) as exception_info:
        OutputSignals([signal1, signal2], {1: 0, 1: 1})  # Dictionary with duplicate keys
    assert str(exception_info.value) == 'Signal UIDs must in fact be unique'

# Test 4 - passed
def test_get_signal_from_uid(sample_signals, sample_graph):
    """Test get_signal_from_uid method."""
    signal_in, signal_1, signal_2, signal_3 = sample_signals
    out_signals = OutputSignals([signal_in, signal_1, signal_2, signal_3], sample_graph)  # added signal_3
    assert out_signals.get_signal_from_uid(1) == signal_1
    assert out_signals.get_signal_from_uid(2) == signal_2
    assert out_signals.get_signal_from_uid(3) == signal_3

# Test 5 - passed
def test_get_source_uid(sample_signals, sample_graph):
    """Test get_source_uid method."""
    signal_in, signal_1, signal_2, signal_3 = sample_signals
    out_signals = OutputSignals([signal_in, signal_1, signal_2, signal_3], sample_graph)
    assert out_signals.get_source_uid(1) == -1  # Parent of signal_1 (UID 1) is signal_in (graph ID 0)
    assert out_signals.get_source_uid(2) == -1  # Parent of signal_2 (UID 2) is signal_in (graph ID 0)
    assert out_signals.get_source_uid(3) == 2  # Parent of signal_3 (UID 3) is signal_2 (graph ID 2)

# Test 6 - passed
def test_outputSignals_properties(sample_signals, sample_graph):
    """Test OutputSignals properties."""
    signal_in, signal_1, signal_2, signal_3 = sample_signals
    out_signals = OutputSignals([signal_in, signal_1, signal_2, signal_3], sample_graph)
    assert out_signals.insignal == signal_in
    assert out_signals.signal_length == len(signal_in)
    assert out_signals.signal_fs == signal_in.fs
    assert out_signals.outsignals == [signal_1, signal_2, signal_3]
    assert out_signals.graph == sample_graph
    assert len(out_signals) == 3

# Test 7 - passed
def test_aLL_plots_with_empty_outSignals( sample_signals, sample_graph): # note: the input signal is the first signal in the list, and the rest of the signals are the output signals.
    """Test correlogram, autocorrelates and correlate_with with empty outsignals."""
    signal_in, signal_1, signal_2, signal_3 = sample_signals
    out_signals = OutputSignals([signal_in], sample_graph)  # Only input signal
    # correlogram
    with pytest.raises(ValueError) as exception_info:
        out_signals.correlogram()
    assert str(exception_info.value) == "Number of rows must be a positive integer, not -1"
    # autocorrelates
    with pytest.raises(ValueError) as exception_info:
        out_signals.autocorrelates()
    assert str(exception_info.value) == 'Number of rows must be a positive integer, not 0'
    # correlate_with
    with pytest.raises(ValueError) as exception_info:
        out_signals.correlate_with(signal_in)
    assert str(exception_info.value) == 'Number of rows must be a positive integer, not 0'

# Test 8 - passed
def test_correlogram_single_output(sample_signals, sample_graph):
    """Test correlogram with a single output signal."""
    signal_in, signal_1, signal_2, signal_3 = sample_signals
    out_signals = OutputSignals([signal_in, signal_1], sample_graph)
    with pytest.raises(ValueError) as e_info:
        out_signals.correlogram()
    assert str(e_info.value) == "Number of rows must be a positive integer, not 0"

# Test 9 - failed (complicated error)
def test_autocorrelates_single_output( sample_signals, sample_graph):
    """Test autocorrelates with a single output signal."""
    signal_in, signal_1, signal_2, signal_3 = sample_signals
    out_signals = OutputSignals([signal_in, signal_1], sample_graph)
    out_signals.autocorrelates()

# Test 10 - failed (complicated error)
def test_correlate_with_single_output(sample_signals, sample_graph):
    """Test correlate_with with a single output signal."""
    signal_in, signal_1, signal_2, signal_3 = sample_signals
    out_signals = OutputSignals([signal_in, signal_1], sample_graph)
    out_signals.correlate_with(signal_in)