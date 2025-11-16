class RootedTree:
  def __init__(self): # root id = 0; node ids = 1, 2, 3...
    self.num_nonroot = 0
    self.parent = [0]
    self.child = [[]]

  def add(self, source_id=0):
    if source_id not in range(self.num_nonroot+1): raise Exception()
    self.num_nonroot += 1
    node_id = self.num_nonroot
    self.parent += [source_id]
    self.child += [[]]
    self.child[source_id] += [node_id]

  def add_parallel(self, n, source_id=0):
    if source_id not in range(self.num_nonroot+1): raise Exception()
    for _ in range(n):
      self.add(source_id=source_id)

  def is_parallel(self):
    return all(self.parent[i]==0 for i in range(1, self.num_nonroot+1))

  def add_series(self, n, source_id=0):
    if source_id not in range(self.num_nonroot+1): raise Exception()
    parent = source_id
    for _ in range(n):
      self.add(source_id=parent)
      parent = self.num_nonroot

  def is_series(self):
    return all(self.parent[i]==i-1 for i in range(1, self.num_nonroot+1))

  def propagate_down(self, root_output, func): # func(node_input, node_id) -> node_output
    outputs = [None for _ in range(self.num_nonroot+1)]
    outputs[0] = root_output
    q = sum(self.child, [])
    for p in q:
      outputs[p] = func(outputs[self.parent[p]], p)
    return outputs