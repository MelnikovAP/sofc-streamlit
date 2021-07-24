import numpy as np
import streamlit as st
import requests
import pandas as pd
from matplotlib import pyplot as plt


API_URL = 'http://178.154.215.108/tote/ivc'
works = ['Introduction', '1. Getting IVC', '2. The second work']

def get_ivc(payload):
    result = requests.post(API_URL, data=payload)
    return result.json()

# ==============================
# Sidebar
# ==============================

tab_selected = st.sidebar.selectbox('Select the work:', works)

if tab_selected != 'Introduction':
    st.sidebar.subheader('Model parameters')
    temperature = st.sidebar.slider(
        'Temperature, '+u'\N{DEGREE SIGN}'+'C', 
        750, 1200, 900, step=1
    )
    sigma = st.sidebar.slider(
        u'\u03C3'+', sigma', 1.8, 2.2, 2.0, step=.01
    )
    ethick = st.sidebar.slider(
        'Electrolyte thickness, µm', 15, 55, 50, step=1
    )
    if st.sidebar.button('Calculate'):
        payload = {
            'temperature': temperature,
            'sigma': sigma,
            'ethick': ethick * pow(10, -6)
        }
        result = get_ivc(payload)
        fig, ax = plt.subplots()
        ax.plot(result['i'], result['Eload'])
        st.pyplot(fig)
# ==============================

# ==============================
# Main page
# ==============================

if tab_selected == 'Introduction':
    col1, col2 = st.beta_columns(2)
    with col1:
        st.image('./scheme.gif')
    with col2:
        '''
        ![](https://raw.githubusercontent.com/MelnikovAP/sofc-streamlit/master/Images/oxygen.png) hydrogen  
        ![](https://raw.githubusercontent.com/MelnikovAP/sofc-streamlit/master/Images/hydrogen.png) oxygen
        ![](https://raw.githubusercontent.com/MelnikovAP/sofc-streamlit/master/Images/water.png)
        ![](https://raw.githubusercontent.com/MelnikovAP/sofc-streamlit/master/Images/electron.png)
        '''




if 'journal' not in st.session_state:
    st.session_state.journal = []

# 2 --- you can add some css to your Streamlit app to customize it
# TODO: Change values below and observer the changes in your app
# st.markdown(
#         f"""
# <style>
#     .reportview-container .main .block-container{{
#         max-width: 90%;
#         padding-top: 5rem;
#         padding-right: 5rem;
#         padding-left: 5rem;
#         padding-bottom: 5rem;
#     }}
#     img{{
#     	max-width:40%;
#     	margin-bottom:40px;
#     }}
# </style>
# """,
#         unsafe_allow_html=True,
#     )
#######################################





# 3 --- build the structure of your app


# Streamlit apps can be split into sections


# container -> horizontal sections
# columns -> vertical sections (can be created inside containers or directly in the app)
# sidebar -> a vertical bar on the side of your app


# here is how to create containers
header_container = st.beta_container()
stats_container = st.beta_container()	
#######################################



# You can place things (titles, images, text, plots, dataframes, columns etc.) inside a container
#with header_container:
