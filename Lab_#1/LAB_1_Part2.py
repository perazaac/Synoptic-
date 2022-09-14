#!/usr/bin/env python
# coding: utf-8

# ## Lab 1 Part II: Accessing Rawinsonde Data
# #### 9/7/2022
# 
# 
# The following tutorial is the second of three parts of the python portion of Lab 1.  In this part we will focus on how to work with rawinsonde data in python using the modules Scipy and Pandas.  
# <br />
# 
# ### Module Documentation
# 1. Pandas DataFrame: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html
# 2. The datetime function from the datetime module: https://docs.python.org/3/library/datetime.html
# 3. Wyoming Upper Air Data via Siphon: https://unidata.github.io/siphon/latest/api/simplewebservice.html#module-siphon.simplewebservice.wyoming
# 
# 
# <br /><br />
# 
# If you have any questions about the code below always feel free to reach out to me at mpvossen@uwm.edu and I am always willing to further explain the code. <br /> <br />
# 
# 
# 1.  Before we start downloading rawinsonde data we first need to import the python modules that we need
# <br />

# In[1]:


#from python's date and time module (from datetime) import the ability to work with date and times (import datetime)
from datetime import datetime

#using the module siphon and its ability to retrieve files from online (.simplewebservice) specifically for the university of wyoming (.wyoming), 
#import the ability to download the university of wyoming's upper air data.
from siphon.simplewebservice.wyoming import WyomingUpperAir


# <br /><br />
# 2. Now with the packages we need, lets choose what rawinsonde observation we want.  Below I set the date, time, and location of the rawinsonde observation that I want.  For the assignment we will want Green Bay's (GRB) observation for 9/7/2022.

# In[2]:


#Here I set the create the variable that holds the time in UTC that we want the rawinsonde data for.
#datetime(year, month, day, hour)
sounding_date = datetime(2022,9,7,12)

#We want to set a variable to the name of the rawinsonde station we want the observation for.
station = "GRB"


# <br /><br />
# 3. With our location and time of the observation set we are ready to download the rawinsonde file using siphon.

# In[ ]:


#using the wyoming upper air rawinsonde downloader we imported above (WyomingUpperAir), get the data (.request_data) for the location and time we specifed before.  
#also with the .set_index("pressure") I set the index of the data frame to be pressure so we can use data.loc[pressure] to get the data at a specified pressure.
upper_air_data = WyomingUpperAir.request_data(sounding_date, station).set_index("pressure")

#display the rawinsonde observation.  Once again the data is in a pandas dataframe
upper_air_data


# <br /><br />
# 3. To complete the next part we need to see what the pandas dataframe column names are.  By looking at the column names we are also looking at what observation variables are contained in the rawinsonde data

# In[ ]:


#display the column names for the rawinsonde data dataframe.
upper_air_data.columns


# <br /><br />
# 4. Just like with the METAR data, we are able to parse out specific variables for specific heights using the pandas dataframe syntax.

# In[ ]:


#display the geopotential height at 500 hPa for the rawinsonde observation we downloaded
upper_air_data["height"].loc[500]


# <br /><br />
# 5.  In the code section below, parse out the temperature, dewpoint, wind speed, wind direction, and geopotential height for 300 hPa.  Display all of the variables so you can use them to answer question 8 in lab 1.

# In[ ]:


from datetime import datetime
from siphon.simplewebservice.wyoming import WyomingUpperAir
sounding_date = datetime(2022,9,7,12)
station = "GRB"
upper_air_data = WyomingUpperAir.request_data(sounding_date, station).set_index("pressure")
upper_air_data
upper_air_data.columns
upper_air_data["height"].loc[300]
upper_air_data["temperature"].loc[300]
upper_air_data["dewpoint"].loc[300]
upper_air_data["speed"].loc[300]
upper_air_data["direction"].loc[300]


# <br /><br />
# 
# You have now completed part II of the python portion of the lab.  Be sure to submit the fully rendered Jupyter Notebook on GitHub when you are finished.
# 
