import pytest
from src.GEFs_core.Filter import Filter
from src.GEFs_core.Signal import Signal
from src.GEFs_core.FilterBank import FilterBank

#======================
# Filter Fixture
#======================
# (can be extended to include more filters)
@pytest.fixture
def sample_filter():
    """Provides a sample Filter object."""
    return Filter(Ap=0.1, bp=1, Bu=3)

#======================
# Signal Fixture
#======================
# (can be extended to include more signals)
@pytest.fixture
def sample_signal():
    """Provides a sample Signal object."""
    return Signal(data=[1, 2, 3], fs=1)

#======================
# Filter Bank Fixtures
#======================
# (can be extended to include more filter banks)
# Parallel
@pytest.fixture
def sample_parallel_filter_bank():
    """Provides a FilterBank with parallel filters."""
    return FilterBank(topology='parallel', Ap=[0.1, 0.2], bp=[1, 2], Bu=[3, 4])

# Series
@pytest.fixture
def sample_series_filter_bank():
    """Provides a FilterBank with series filters."""
    return FilterBank(topology='series', Ap=[0.1, 0.2], bp=[1, 2], Bu=[3, 4])

#======================
# Core Tests
#======================

# Test 1 - Passed
def test_filterbank_initialisation_parameter_combinations():
    """Tests FilterBank initialization with various parameter combinations and checks the number of filters created accordingly."""
    test_data = [
        ({'Ap': 0.2, 'bp': 2, 'Bu': 4}, 1), #length of variable
        ({'topology': 'parallel', 'Ap': [0.3, 0.4], 'bp': [3, 4], 'Bu': [5, 6]}, 2),
        ({'topology': 'series', 'Ap': [0.5, 0.6], 'bp': 5, 'Bu': [7, 8]}, 2)
    ]
    for params, expected_count in test_data:
        fb = FilterBank(**params)
        print(f"Testing with parameters: {params}")
        print(f"Expected filter count: {expected_count}")
        print(f"Actual filter count: {len(fb.filters)}")
        assert len(fb.filters) == expected_count
        print("Test passed.\n")

# Test 2 - passed
def test_filterbank_initialization_preexisting_filters(sample_filter):
    """Tests FilterBank initialization with pre-existing Filter objects."""
    f2 = Filter(Ap=0.2, bp=2, Bu=4)
    fb = FilterBank(filters=[sample_filter, f2], topology='parallel')
    assert len(fb.filters) == 2

# Test 3 - failed (why? it assume that there is 3 filters inside the bank)
def test_filterbank_graph_parallel_topology_parents(sample_parallel_filter_bank):
    """Tests the parent structure of a parallel FilterBank graph."""
    assert sample_parallel_filter_bank.graph.parent == [0, 0]

# Test 4 - failed (why? it assume that there is 3 filters inside the bank)
def test_filterbank_graph_series_topology_parents(sample_series_filter_bank):
    """Tests the parent structure of a series FilterBank graph."""
    assert sample_series_filter_bank.graph.parent == [0, 1]

# Test 5 - passed
def test_filterbank_initialization_invalid_topology_raises_error():
    """Tests that an invalid topology raises an error during initialization."""
    with pytest.raises(Exception, match='Invalid topology name'):
        FilterBank(topology='invalid', Ap=[0.1, 0.2], bp=[1, 2],Bu=4)

# Test 6 - failed (it assume its 3 filters while its 2)
def test_filterbank_graph_add_filter_to_empty_bank(sample_filter):
    """Tests adding a filter to an empty FilterBank."""
    fb = FilterBank(topology='parallel', Ap=[0.1, 0.2], bp=[1, 2],Bu=4)
    fb.add(sample_filter)
    assert len(fb.filters) == 2

# Test 7 - failed ( needs an exception handling)
def test_filterbank_graph_add_filter_with_invalid_source_raises_error(sample_filter):
    """Tests adding a filter with an invalid source UID raises an error."""
    fb = FilterBank(topology='parallel', Ap=[0.1, 0.2], bp=[1, 2], Bu=4)
    fb.add(sample_filter, source_uid=999) # cause an error

# Test 8 - passed
def test_filterbank_signal_processing_parallel(sample_parallel_filter_bank, sample_signal):
    """Tests signal processing through a parallel FilterBank."""
    output = sample_parallel_filter_bank.process_signal(sample_signal)
    assert len(output.outsignals) == 2
    assert all(isinstance(sig, Signal) for sig in output.outsignals)

# Test 9 - passed
def test_filterbank_signal_processing_series(sample_series_filter_bank, sample_signal):
    """Tests signal processing through a series FilterBank."""
    output = sample_series_filter_bank.process_signal(sample_signal)
    assert len(output.outsignals) == 2
    assert output.outsignals[1].mode_t is not None

# Test 10 - failed (2 errors here needs exception handling)
def test_filterbank_signal_processing_empty_bank_raises_error(sample_signal):
    """Tests that signal processing through an empty FilterBank raises an error."""
    filter_bank = FilterBank() # could not handle empty initialization
    filter_bank.process_signal(sample_signal) # cause an error that expose a lot of system details

# Test 11 - passed
def test_filterbank_error_missing_cf_for_unnormalized_raises_error():
    """Tests that missing 'cf' for unnormalized filters raises an error."""
    with pytest.raises(Exception, match='Please provide characteristic frequency of filter'):
        FilterBank(topology='parallel', fpeak=[10], Nf=[10])

# Test 12 - passed
def test_filterbank_relationships_source_identification():
    """Tests the identification of source filters in a FilterBank."""
    f1 = Filter(Ap=0.1, bp=1, Bu=1)
    f2 = Filter(Ap=0.2, bp=2, Bu=2)
    f1.uid = 1
    f2.uid = 2
    fb = FilterBank(filters=[f1, f2], topology='series')
    assert fb.get_source_uid(2) == 1

# Test 13 - passed
def test_filterbank_relationships_child_filters():
    """Tests the identification of child filters in a FilterBank."""
    f1 = Filter(Ap=0.1, bp=1, Bu=1)
    f2 = Filter(Ap=0.2, bp=2, Bu=2)
    f1.uid = 1
    f2.uid = 2
    fb = FilterBank(filters=[f1, f2], topology='series')
    assert fb.get_uids_fed_into(1) == [2]

# Test 14 - passed
def test_filterbank_relationships_uid_mapping():
    """Tests the UID mapping within a FilterBank."""
    f1 = Filter(Ap=0.1, bp=1, Bu=1)
    f2 = Filter(Ap=0.2, bp=2, Bu=2)
    f1.uid = 1
    f2.uid = 2
    fb = FilterBank(filters=[f1, f2], topology='series')
    assert fb._uid2graphid[1] == 1
    assert fb._graphid2uid[1] == 1