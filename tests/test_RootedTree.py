import pytest
from src.GEFs_core.RootedTree import RootedTree

#======================
# Code Tests
#======================

# Test 1 - passed
def test_rootedtree_initialization():
    """Tests RootedTree initialization."""
    tree = RootedTree()
    assert tree.num_nonroot == 0 # number of non root nodes
    assert tree.parent == [0] # The root node (node 0) has itself as its parent
    assert tree.child == [[]]

# Test 2 - passed
def test_add_node():
    """Tests adding a single node to the tree."""
    tree = RootedTree()
    tree.add()
    assert tree.num_nonroot == 1
    assert tree.parent == [0, 0] # parent attribute after adding a node
    assert tree.child == [[1], []] # list of lists, where tree.child[i] is a list of the children of node

# Test 3 - passed
def test_add_node_invalid_source_id():
    """Tests adding a node with an invalid source ID."""
    tree = RootedTree()
    with pytest.raises(Exception):
        tree.add(source_id=1)  #not in range self.num_nonroot+1 (same happening with add_parallel and add_series)

# Test 4 - passed
def test_add_zero_nodes():
    """Tests the behavior of RootedTree class when adding zero nodes in parallel and series."""
    tree = RootedTree()
    tree.add_parallel(0)
    tree.add_series(0)
    assert tree.num_nonroot == 0

# Test 5 - passed
def test_add_parallel():
    """Tests adding multiple nodes in parallel."""
    tree = RootedTree()
    tree.add_parallel(3)
    assert tree.num_nonroot == 3
    assert tree.parent == [0, 0, 0, 0]
    assert tree.child == [[1, 2, 3], [], [], []]

# Test 6 - passed
def test_add_series():
    """Tests adding multiple nodes in series."""
    tree = RootedTree()
    tree.add_series(3)
    assert tree.num_nonroot == 3
    assert tree.parent == [0, 0, 1, 2]
    assert tree.child == [[1], [2], [3], []]

# Test 7 - passed
def test_is_parallel_true():
    """Tests if the tree is parallel when it should be."""
    tree = RootedTree()
    tree.add_parallel(3)
    assert tree.is_parallel() is True

# Test 8 - passed
def test_is_parallel_false():
    """Tests if the tree is not parallel when it shouldn't be."""
    tree = RootedTree()
    tree.add_series(3)
    assert tree.is_parallel() is False

# Test 9 - passed
def test_is_series_true():
    """Tests if the tree is series when it should be."""
    tree = RootedTree()
    tree.add_series(3)
    assert tree.is_series() is True

# Test 10 - passed
def test_is_series_false():
    """Tests if the tree is not series when it shouldn't be."""
    tree = RootedTree()
    tree.add_parallel(3)
    assert tree.is_series() is False

# Test 11 - passed
def test_large_scale_operations():
    """ Tests the ability to handle large-scale series operations."""
    tree = RootedTree()
    tree.add_series(100)
    assert tree.is_series() is True
    assert len(tree.parent) == 101

# Test 12 - passed
def test_propagate_down():
    """Tests propagating values down the tree."""
    tree = RootedTree()
    tree.add_parallel(2)
    tree.add_series(2, source_id=1)# add 2 nodes under node 1
    tree_outputs = tree.propagate_down(
        root_output=10,
        func=lambda node_input, node_id: node_input + node_id
    )
    assert tree_outputs == [10, 11, 12, 14, 18]

# Test 13 - passed
def test_propagate_order():
    """ Tests the propagate_down method to ensure it traverses the tree in the correct order."""
    tree = RootedTree()
    tree.add_parallel(3)
    outputs = tree.propagate_down(0, lambda x, nid: f"{x}{nid}")
    assert outputs == [0, '01', '02', '03'] # visiting the first child of the root. The lambda function concatenates '0' (the accumulated value) with '1' (the node ID).