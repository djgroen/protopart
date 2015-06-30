import numpy as np

class CSR:
  """ Class definition for the CSR format, which is directly adopted by ParMETIS and PPStee.
  """
  def __init__(self):
    
    # Official CSR fields
    self.vtxdist = []
    self.xadj = []
    self.adjncy = []
    self.wgts = []
    
    # Added fields to aid in analysis and visualization
    self.vids = []
    self.vcoords = []