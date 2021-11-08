#!/usr/bin/env python
# coding: utf-8

# In[1]:


##get_ipython().run_line_magic('load_ext', 'autoreload')
##get_ipython().run_line_magic('autoreload', '2')
##get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


import numpy as np
from conrad import *


# In[4]:


# dimensions
voxels = 1000
beams = 200

# voxel labels
labels = (np.random.rand(voxels) > 0.2).astype(int)

# dose matrix with target voxels receiving ~3x radiation of non-target voxels
# label 0 = TUMOR, label 1 = OAR
A = np.random.rand(voxels, beams)
FACTOR = 3
for i, label in enumerate(labels):
    if label == 0:
        A[i, :] *= FACTOR
        
# construct case
case = Case()
case.anatomy = Anatomy([Structure(0, 'TUMOR', True), Structure(1, 'OAR', False)])
case.physics = Physics(dose_matrix=A, voxel_labels=labels)


# initialize CasePlotter object
graphics = CasePlotter(case)


# In[5]:


# add DVH constraints
case.anatomy['TUMOR'].constraints += D(15) < 1.05 * Gy
case.anatomy['TUMOR'].constraints += D(85) > 0.9 * Gy
case.anatomy['OAR'].constraints += D(50) < 0.4 * Gy


# In[6]:


# solve: WITHOUT slack, ONE pass (restricted DVH constraints)
status, run = case.plan(use_slack=False)
print('SOLVER CONVERGED?', status)
graphics.plot(run.plotting_data, plotfile='DVH_Constraints1.png')


# In[7]:


# solve: WITHOUT slack, TWO passes (exact DVH constraints)
status, run = case.plan(use_slack=False, use_2pass=True)
print('SOLVER CONVERGED?', status)
graphics.plot(run.plotting_data, second_pass=True, plotfile='DVH_Constraints2.png')


# In[8]:


# additional DVH constraints makes no-slack problem infeasible
case.anatomy['TUMOR'].constraints += D(99) > 0.99 * Gy
case.anatomy['TUMOR'].constraints += D(1) < 1.01 * Gy


# In[9]:


# solving without slacks will result in infeasibility
status, run = case.plan(dvh_slack=False)
print('SOLVER CONVERGED?', status)


# In[10]:


# solve: WITH slack, ONE pass
# (N.B.: use_slack=True by default, so could call w/o use_slack argument)
status, run = case.plan(use_slack=True)
print('SOLVER CONVERGED?', status)
graphics.plot(run.plotting_data, plotfile='DVH_Constraints3.png')


# In[11]:


# solve: WITH slack, TWO passes 
status, run = case.plan(use_2pass=True)
print('SOLVER CONVERGED?', status)
graphics.plot(run.plotting_data, second_pass=True, plotfile='DVH_Constraints4.png')


# In[ ]:




