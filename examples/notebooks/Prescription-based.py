#!/usr/bin/env python
# coding: utf-8

# In[ ]:


##get_ipython().run_line_magic('load_ext', 'autoreload')
##get_ipython().run_line_magic('autoreload', '2')
##get_ipython().run_line_magic('matplotlib', 'inline')


# In[ ]:


from numpy.random import rand
from conrad import *


# In[ ]:


# prescription, as list of dictionaries:
rx = [{
        'name': 'tumor',
        'label': 0,
        'is_target': True,
        'dose': '30 Gy',
        'constraints': ['D85 > 0.9 rx', 'D99 < 1.1 rx'],
    },{
        'name': 'OAR',
        'label': 1,
        'is_target': False,
        'constraints': ['mean < 5 Gy', 'D95 < 15 Gy'],
    },{
        'name': 'body',
        'label': 2,
        'is_target': False,
        'constraints': ['mean < 10 Gy'],
    }]

# same prescription, as object
prescription = Prescription(rx)


# In[ ]:


TUMOR = 0
OAR = 1
BODY = 2

# dimensions, voxel labels and dose matrix 
voxels, beams = 2000, 300
v = (rand(voxels) > 0.3).astype(int) + (rand(voxels) > 0.3).astype(int)
A = rand(voxels, beams)
factor = 3.
for voxel, label in enumerate(v):
    if label == TUMOR:
        A[voxel, :] *= factor
        
physics = Physics(dose_matrix=A, voxel_labels=v)


# In[ ]:


# case from prescription (take structures, suppress constraints)
case = Case(physics=physics, prescription=rx, suppress_rx_constraints=True)
graphics = CasePlotter(case)
print(case.anatomy)


# In[ ]:


_, run = case.plan()
graphics.plot(run.plotting_data, plotfile="Prescription-based1.png")


# In[ ]:


# transfer constraints from prescription to anatomy
case.transfer_rx_constraints_to_anatomy()
print(case.anatomy)
_, run = case.plan()
graphics.plot(run.plotting_data, plotfile="Prescription-based2.png")


# In[ ]:


# alternate syntax 1 (don't suppress constraints)
case = Case(physics=physics, prescription=rx)
print(case.anatomy)
_, run = case.plan()
graphics.plot(run.plotting_data, plotfile="Prescription-based3.png")


# In[ ]:


# alternate syntax 2 (initialize with prescription object)
case = Case(physics=physics, prescription=prescription)
print(case.anatomy)
_, run = case.plan()
graphics.plot(run.plotting_data, plotfile="Prescription-based4.png")


# In[ ]:


# alternate syntax 3 (initialize without prescription, add later as list of dictionaries)
case = Case(physics=physics)
case.prescription = rx

# N.B.: automatic transfer of constraints from prescription to anatomy only happens in
# constructor. Explicit call of transfer() method required when adding prescription to 
# an already-initialized case 
print(case.anatomy)
_, run = case.plan()
graphics.plot(run.plotting_data, plotfile="Prescription-based5.png")


# In[ ]:


# transfer constraints and re-plan
case.transfer_rx_constraints_to_anatomy()
print(case.anatomy)
_, run = case.plan()
graphics.plot(run.plotting_data, plotfile="Prescription-based6.png")


# In[ ]:


# alternate syntax 4 (initialize without prescription, add later as object)
case = Case(physics=physics)
case.prescription = prescription
print(case.anatomy)
_, run = case.plan()
graphics.plot(run.plotting_data, plotfile="Prescription-based7.png")


# In[ ]:


# transfer constraints and re-plan
case.transfer_rx_constraints_to_anatomy()
print(case.anatomy)
_, run = case.plan()
graphics.plot(run.plotting_data, plotfile="Prescription-based8.png")

