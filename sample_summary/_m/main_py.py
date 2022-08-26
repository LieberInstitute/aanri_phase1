#!/usr/bin/env python
# coding: utf-8

# # Pie chart of samples

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


df = pd.read_csv("../../input/phenotypes/merged/_m/merged_phenotypes.csv", index_col=0)
df = df[df["Race"].isin(["AA", "CAUC"]) & (df["Dx"].isin(["Control", "Schizo"]))].copy()
df.Race = df.Race.astype("category").cat.rename_categories({'CAUC': 'EA'})
df.head()


# In[4]:


df.groupby("Region").size()


# In[5]:


data = df.groupby(["Region", "Race"]).size().reset_index().rename(columns={0:"N"})
data


# In[6]:


caudate = data[(data["Region"] == "Caudate")].drop("Region", axis=1).set_index("Race")
dlpfc = data[(data["Region"] == "DLPFC")].drop("Region", axis=1).set_index("Race")
gyrus = data[(data["Region"] == "DentateGyrus")].drop("Region", axis=1).set_index("Race")
hippocampus = data[(data["Region"] == "HIPPO")].drop("Region", axis=1).set_index("Race")


# In[7]:


caudate.N.plot.pie(autopct="%.1f%%")
plt.savefig('caudate_pie.pdf')


# In[8]:


dlpfc.N.plot.pie(autopct="%.1f%%");
plt.savefig('dlpfc_pie.pdf')


# In[9]:


gyrus.N.plot.pie(autopct="%.1f%%");
plt.savefig('dentate_gyrus_pie.pdf')


# In[10]:


hippocampus.N.plot.pie(autopct="%.1f%%")
plt.savefig('hippocampus_pie.pdf')


# In[ ]:




