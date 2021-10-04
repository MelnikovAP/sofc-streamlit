import os
import sys
from pathlib import WindowsPath
import numpy as np
from requests.api import head
import streamlit as st
import requests
import pandas as pd
from bokeh.models import ColumnDataSource, Label, LabelSet, Range1d
from bokeh.plotting import figure, output_file, show
import time

import matplotlib.pyplot as plt

from PIL import Image

# Function to Read and Manupilate Images
def load_image(img):
    im = Image.open(img)
    image = np.array(im)
    return image

# sys.path.append(os.path.dirname(__file__))
# def local_css(file_name):
#     with open(file_name) as f:
#         st.markdown('<style>{}</style>'.format(f.read()), unsafe_allow_html=True)

# local_css(os.path.join(os.path.dirname(__file__),"style.css"))



#st.set_page_config(page_title=None, page_icon=None, layout='wide', initial_sidebar_state='auto')

API_URL = 'http://178.154.215.108/tote/ivc'
# works_eng = ['Introduction',
#          'Work #1: Understanding IVC',
#          'Work #2: How fuel cell works: Thermodynamics',
#          'Work #3: How fuel cell works: Ohmic losses',
#          'Work #4: How fuel cell works: Mass losses',
#          'Work #5: From cell to stack',
#          'Futher investigation']

# works_ru = ['Введение',
#          'Work #1: Understanding IVC',
#          'Work #2: How fuel cell works: Thermodynamics',
#          'Work #3: How fuel cell works: Ohmic losses',
#          'Work #4: How fuel cell works: Mass losses',
#          'Work #5: From cell to stack',
#          'Futher investigation']

# ==============================
# Sidebar
# ==============================

# st.sidebar.image('./Images/logo.png')
# st.sidebar.image(
#     'https://raw.githubusercontent.com/MelnikovAP/pv-streamlit/main/Images/logo.png')




#img_logo = '''
#<img alt='logo' src='https://raw.githubusercontent.com/MelnikovAP/pv-streamlit/main/Images/logo.png' width="100%">
#'''


# img_logo = "<img alt='logo' src='"+ os.path.join(os.path.dirname(__file__),"Images","logo.png") + "' width='100%'>"
# st.sidebar.markdown(img_logo, unsafe_allow_html=True)


hide_full_screen = '''
<style>
[data-testid='stSidebar'] button[title='View fullscreen']  {visibility: hidden;}
</style>
'''
st.sidebar.markdown(hide_full_screen, unsafe_allow_html=True) 

st.sidebar.image(load_image(os.path.dirname(__file__) + r"\Images\logo.png"))





langs = ["RU","ENG"]
col1, col2 = st.sidebar.columns([2,1])
col1.markdown(r"<div style='margin-top:55px;text-align:right;'><hr /></div>", unsafe_allow_html=True)
lang_selected = col2.selectbox("", langs)
#col1.markdown('***')




works = []
# if lang_selected == langs[0]:
#     works = works_eng
# if lang_selected == langs[1]:
#     works = works_ru  

td = {
    "ENG" :[
        [
            'Introduction',
#            'Understanding IVC',
            'Work #1: How SOFC works: Thermodynamics',
            'Work #2: How SOFC works: Ohmic losses',
            'Work #3: How SOFC works: Mass losses',
            'Work #4: From cell to stack',
#            'Futher investigation'
            ],
        'SOFC cloud model',
        'Select the work:',
        'Model parameters',
        'temperature',
        'pressure',
        'electrolyte thickness',
        'conductivity',
        'atm',
        'mA/cm²',
        'S/m',
        'µm',
        'Reset',
        'Run'+u'\u00a0'+'simulation',
        'Allow to refresh',
        'V',
        'mW/cm²',
        'Voltage',
        'Current density',
        'Power',
        "Total consumed H₂, [kg]",
        "Efficiency factor, [%]",
        "Generated power, [kWh]",
        "Generated H₂O, [kg]",
    ],
    "RU" :[
        [
            'Введение',
#            'Understanding IVC',
            'Работа #1: Термодинамика процессов',
            'Работа #2: Омические потери',
            'Работа #3: Концентрационные потери',
            'Работа #4: Калькулятор стека',
#            'Futher investigation'
            ],
        'Модель ТОТЭ',
        'Выбор работы:',
        'Параметры модели',
        'температура',
        'давление',
        'толщина электролита',
        'удельная проводимость',
        'атм',
        'мА/см²',
        'См/м',
        'мкм',
        'Сброс',
        'Запуск'+u'\u00a0'+'симуляции',
        'Обновлять график',
        'В',
        'мВт/см²',
        'Напряжение',
        'Плотность тока',
        'Мощность',
        "Необходимое количество H₂, [кг]",
        "КПД, [%]",
        "Электрическая мощность установки, [кВтч]",
        "Количество H₂O, [кг]",
    ],
}


works = td[lang_selected][0]

st.sidebar.header(td[lang_selected][1])

tab_selected = st.sidebar.selectbox(td[lang_selected][2], works)


temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc = 1000, 1.1, 2.0, 50, 80, 0.95, 0.05, 0.95




if tab_selected in works[1:4]:
    st.sidebar.subheader(td[lang_selected][3])
    temperature = st.sidebar.slider(
        #'T [temperature], '+u'\N{DEGREE SIGN}'+'K',
        'T [' + td[lang_selected][4] + '], K',
        700, 1500, 1000, step=10, key="sld_temperature",
    )

    pressure = st.sidebar.slider(
        'p [' + td[lang_selected][5] + '], ' + td[lang_selected][8], 1.0, 3.0, 1.1, step=.1, key="sld_pressure"
    )

    sigma = st.sidebar.slider(
        u'\u03C3'+' [' + td[lang_selected][7] + '], ' + td[lang_selected][10], 0.0, 3.0, 2.0, step=.2, key="sld_sigma"
    )
    ethick = st.sidebar.slider(
        'd [' + td[lang_selected][6] + '], ' + td[lang_selected][11], 15, 155, 50, step=5, key="sld_ethick"
    )


    jm = st.sidebar.slider(
        'jₘ , ' + td[lang_selected][9], 10, 150, 80, step=5, key="sld_jm"
    )

    H2ac = st.sidebar.slider('H₂, %', 0.8, 0.99,
                             0.97, step=0.01, key="sld_H2ac")
    H2Oac = st.sidebar.slider(
        'H₂O, %', 0.01, 0.07, 0.03, step=0.01, key="sld_H2Oac")


    # O2cc = st.sidebar.slider('O₂, %', 0.01, 1.0, 0.95,
    #                           step=0.01, key="sld_O2cc")


    if "sld_temperature" not in st.session_state:
        st.session_state["sld_temperature"] = 1000
    if "sld_pressure" not in st.session_state:
        st.session_state["sld_pressure"] = 1.1
    if "sld_sigma" not in st.session_state:
        st.session_state["sld_sigma"] = 2.0
    if "sld_ethick" not in st.session_state:
        st.session_state["sld_ethick"] = 50
    if "sld_jm" not in st.session_state:
        st.session_state["sld_jm"] = 80
    if "sld_H2ac" not in st.session_state:
        st.session_state["sld_H2ac"] = 0.97
    if "sld_H2Oac" not in st.session_state:
        st.session_state["sld_H2Oac"] = 0.03
    if "sld_O2cc" not in st.session_state:
        st.session_state["sld_O2cc"] = 0.95

    def _update_sliders(temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc):
        st.session_state["sld_temperature"] = temperature
        st.session_state["sld_pressure"] = pressure
        st.session_state["sld_sigma"] = sigma
        st.session_state["sld_ethick"] = ethick
        st.session_state["sld_jm"] = jm
        st.session_state["sld_H2ac"] = H2ac
        st.session_state["sld_H2Oac"] = H2Oac
        st.session_state["sld_O2cc"] = O2cc

    col1, col2, col3 = st.sidebar.columns([1, 2, 1])
    btn_reset = col1.button(td[lang_selected][12], on_click=_update_sliders, kwargs={
                            "temperature": 1000, "pressure": 1.1, "sigma": 2.0, "ethick": 50, "jm": 80, "H2ac": 0.97, "H2Oac": 0.03, "O2cc": 0.95})
    btn_runsimulation = col2.button(td[lang_selected][13])



# ==============================

# def get_ivc(temperature,sigma,ethick):
#     payload = {
#         'temperature': temperature,
#         'sigma': sigma,
#         'ethick': ethick * pow(10, -6)
#     }
#     result = requests.post(API_URL, data=payload)
#     return

# sys.path.append(os.path.dirname(__file__) + "/../../rev_5")


# from tsofc import TSOFC
# 
# def get_ivc(temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc):
#     payload = {
#         "T": temperature,
#         "P": pressure * 101325.0,
#         "sigma": sigma,
#         "ethick": ethick * pow(10, -6),
#         "jm": jm,
#         "anode_gas_X": "H2:{0:.2f}, H2O:{1:.2f}".format(H2ac, H2Oac),
#         "cathode_gas_X": "O2:{0:.2f}, H2O:0.001".format(O2cc)
#     }
#     # print(payload)
#     result = TSOFC(**payload).getPolarizationCurve()
#     # print(result)
#     return {
#         'i': result[:, 0],
#         'E1': result[:, 1],
#         'E2': result[:, 2],
#         'eta_ohmic': result[:, 3],
#         'Eload': result[:, 4],
#         'eta_con': result[:, 5],
#     }


def get_ivc(temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc):
    payload = {
        "temperature": temperature,
        "pressure": pressure * 101325.0,
        "sigma": sigma,
        "ethick": ethick * pow(10, -6),
        "jm": jm,
        "H2ac": H2ac,
        "H2Oac":H2Oac,
        "O2cc":O2cc
    }
    result = requests.post(API_URL, data=payload).json()
    return {k: np.array(result[k]) for k in result.keys()}


# print(get_ivc(temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc))

# def get_ivc(temperature, sigma, ethick, m, n, H2ac, H2Oac, O2cc):
#     payload = {
#         "T": temperature,
#         "sigma": sigma,
#         "ethick": ethick * pow(10, -6),
#         "m": m * pow(10, -5),
#         "n": n * pow(10, -2),
#         "H2ac": H2ac,
#         "H2Oac": H2Oac,
#         "O2cc": O2cc,
#     }
#     result = requests.post(API_URL, data=payload)
#     print(result,payload)
#     return result




# ==============================
# Main page
# ==============================


header_container = st.container()
for item in works:
    if tab_selected == item:
        with header_container:
            st.title(item)
            '''***'''

if tab_selected == works[0] and lang_selected == langs[1]:

    cflag = False

    if cflag:

        st.markdown('''
        Under construction...
        ''')
    else:

        introduction_container = st.container()
        scheme_container = st.container()
        basics_container = st.container()
        theory_container = st.container()
        whatnext_container = st.container()

        with introduction_container:
            '''
            ### Solid oxide fuel cell (SOFC)
            is an electrochemical conversion device 
            that produces electricity directly from oxidizing a fuel. Fuel cells are 
            characterized by their electrolyte material; the SOFC has a solid oxide 
            or ceramic electrolyte.  
            '''
            '''
            Advantages of this class of fuel cells include high combined heat and 
            power efficiency, long-term stability, fuel flexibility, low emissions, 
            and relatively low cost. The largest disadvantage is the high operating 
            temperature which results in longer start-up times and mechanical and 
            chemical compatibility issues.  
            '''

        with scheme_container:
            st.markdown('<p style="text-align: center; padding-right:20%; padding-left:20%"><img src="https://raw.githubusercontent.com/MelnikovAP/sofc-streamlit/master/Images/scheme.gif", width=80%/></p>',
                        unsafe_allow_html=True)
            '''
            >![](https://raw.githubusercontent.com/MelnikovAP/sofc-streamlit/master/Images/oxygen.png) pure $H_2$ or $H_2+CO$   
            >![](https://raw.githubusercontent.com/MelnikovAP/sofc-streamlit/master/Images/hydrogen.png) $O_2$ from air  
            >![](https://raw.githubusercontent.com/MelnikovAP/sofc-streamlit/master/Images/water.png) $H_2O$ with excess fuel  
            >![](https://raw.githubusercontent.com/MelnikovAP/sofc-streamlit/master/Images/electron.png) $e^{-1}$ that forms electric current  
            '''

        with basics_container:
            '''
            SOFCs are a class of fuel cells characterized by the use of a solid 
            oxide material as the electrolyte. SOFCs use a solid oxide electrolyte 
            to conduct negative oxygen ions from the cathode to the anode. The 
            electrochemical oxidation of the hydrogen, carbon monoxide or other 
            organic intermediates by oxygen ions thus occurs on the anode side. 
            More recently, proton-conducting SOFCs (PC-SOFC) are being developed 
            which transport protons instead of oxygen ions through the electrolyte 
            with the advantage of being able to be run at lower temperatures than 
            traditional SOFCs.  
            '''
            '''
            The electrochemical reactions, that are taking place at the SOFC electrodes, 
            are as follows:
            '''
            st.latex(
                r''' \text{cathode}: \frac{1}{2}O_2+2e^- \rightarrow O^{2-}''')
            with st.expander('Anode description'):
                '''
                The ceramic anode layer must be very porous to allow the fuel to flow 
                towards the electrolyte. Consequently, granular matter is often selected 
                for anode fabrication procedures. Like the cathode, it must conduct 
                electrons, with ionic conductivity a definite asset. The anode is 
                commonly the thickest and strongest layer in each individual cell, 
                because it has the smallest polarization losses, and is often the layer 
                that provides the mechanical support. Electrochemically speaking, the 
                anode’s job is to use the oxygen ions that diffuse through the 
                electrolyte to oxidize the hydrogen fuel. The oxidation reaction 
                between the oxygen ions and the hydrogen produces heat as well as 
                water and electricity. If the fuel is a light hydrocarbon, for example, 
                methane, another function of the anode is to act as a catalyst for steam 
                reforming the fuel into hydrogen. This provides another operational 
                benefit to the fuel cell stack because the reforming reaction is 
                endothermic, which cools the stack internally. The most common material 
                used is a cermet made up of nickel mixed with the ceramic material that 
                is used for the electrolyte in that particular cell, typically YSZ (yttria 
                stabilized zirconia). These nanomaterial-based catalysts, help stop the 
                grain growth of nickel. Larger grains of nickel would reduce the contact 
                area that ions can be conducted through, which would lower the cells 
                efficiency. [Perovskite materials](https://en.wikipedia.org/wiki/Perovskite_(structure)) 
                (mixed ionic/electronic conducting ceramics) have been shown to produce a power density 
                of 0.6 W/cm$^2$ at 0.7 V at 800 °C which is possible because they have the ability 
                to overcome a larger activation energy.
                '''
            st.latex(r''' \text{anode}: H_2+O^{2-} \rightarrow H_2O+2e^-''')
            with st.expander('Cathode description'):
                '''
                Cathode materials must be, at a minimum, electrically conductive. Currently, 
                lanthanum strontium manganite (LSM) is the cathode material of choice for 
                commercial use because of its compatibility with doped zirconia electrolytes. 
                Mechanically, it has a similar coefficient of thermal expansion to YSZ and thus 
                limits stress buildup because of CTE mismatch. Also, LSM has low levels of chemical 
                reactivity with YSZ which extends the lifetime of the materials. Unfortunately, LSM 
                is a poor ionic conductor, and so the electrochemically active reaction is limited 
                to the triple phase boundary (TPB) where the electrolyte, air and electrode meet. 
                LSM works well as a cathode at high temperatures, but its performance quickly falls 
                as the operating temperature is lowered below 800 °C. In order to increase the reaction 
                zone beyond the TPB, a potential cathode material must be able to conduct both electrons 
                and oxygen ions. Composite cathodes consisting of LSM YSZ have been used to increase this 
                triple phase boundary length. Mixed ionic/electronic conducting (MIEC) ceramics, such as 
                perovskite LSCF, are also being researched for use in intermediate temperature SOFCs as 
                they are more active and can make up for the increase in the activation energy of the reaction.
                '''

            '''
            Based on these half-cell reactions, the overall reaction occurring is:
            '''

            st.latex(r''' H_2+\frac{1}{2}O_2 \rightarrow H_2O''')

            with st.expander('Electrolyte description'):
                '''
                The electrolyte is a dense layer of ceramic that conducts oxygen ions. Its electronic 
                conductivity must be kept as low as possible to prevent losses from leakage currents. 
                The high operating temperatures of SOFCs allow the kinetics of oxygen ion transport to 
                be sufficient for good performance. However, as the operating temperature approaches the 
                lower limit for SOFCs at around 600 °C, the electrolyte begins to have large ionic transport 
                resistances and affect the performance. Popular electrolyte materials include yttria-stabilized 
                zirconia (YSZ) (often the 8% form 8YSZ), scandia stabilized zirconia (ScSZ) (usually 9 mol% Sc$_2$O$_3$ – 9ScSZ) 
                and gadolinium doped ceria (GDC). The electrolyte material has crucial influence on the cell performances. 
                Detrimental reactions between YSZ electrolytes and modern cathodes such as lanthanum strontium cobalt 
                ferrite (LSCF) have been found, and can be prevented by thin (<100 nm) ceria diffusion barriers.
                If the conductivity for oxygen ions in SOFC can remain high even at lower temperatures (current 
                target in research ~500 °C), material choices for SOFC will broaden and many existing problems 
                can potentially be solved. 
                '''
            with st.expander('Interconnect'):
                '''
                The interconnect can be either a metallic or ceramic layer that sits between each individual cell. 
                Its purpose is to connect each cell in series, so that the electricity each cell generates can be 
                combined. Because the interconnect is exposed to both the oxidizing and reducing side of the cell at 
                high temperatures, it must be extremely stable. For this reason, ceramics have been more successful 
                in the long term than metals as interconnect materials. However, these ceramic interconnect materials 
                are very expensive as compared to metals. Nickel- and steel-based alloys are becoming more promising 
                as lower temperature (600–800 °C) SOFCs are developed. The material of choice for an interconnect in 
                contact with Y8SZ is a metallic 95Cr-5Fe alloy. Ceramic-metal composites called 'cermet' are also 
                under consideration, as they have demonstrated thermal stability at high temperatures and excellent 
                electrical conductivity.
                '''
            '''
            SOFCs operate at very high temperatures, typically between 500 and 
            1,000 °C. At these temperatures, SOFCs do not require expensive 
            platinum catalyst material, as is currently necessary for lower 
            temperature fuel cells such as [PEMFCs](https://en.wikipedia.org/wiki/Proton-exchange_membrane_fuel_cell), and are not vulnerable to 
            carbon monoxide catalyst poisoning. However, vulnerability to 
            sulfur poisoning has been widely observed and the sulfur must be 
            removed before entering the cell through the use of adsorbent beds 
            or other means.
            ***
            '''

        with theory_container:
            '''
            ### Polarizations
            Polarizations, or overpotentials, are losses in voltage due to imperfections in materials, 
            microstructure, and design of the fuel cell. Polarizations result from ohmic resistance of 
            oxygen ions conducting through the electrolyte ($iR\Omega$), electrochemical activation barriers 
            at the anode and cathode, and finally concentration polarizations due to inability of gases 
            to diffuse at high rates through the porous anode and cathode. 
            The cell voltage can be calculated using the following equation:
            '''
            st.latex(r''' V=E_0-iR_\omega-\eta_{cathode}-\eta_{anode}''')
            '''
            where:  
            >$E_0=$ [Nernst potential](https://en.wikipedia.org/wiki/Reversal_potential) of the reactants  
            >$R=$ [Thévenin equivalent](https://en.wikipedia.org/wiki/Thévenin%27s_theorem) 
            resistance value of the electrically conducting portions of the cell  
            >$\eta_{cathode}=$ polarization losses in the cathode  
            >$\eta_{anode}=$ polarization losses in the anode  
            '''
            '''
            In SOFCs, it is often important to focus on the ohmic and concentration polarizations since 
            high operating temperatures experience little activation polarization. However, as the lower 
            limit of SOFC operating temperature is approached (~600 °C), these polarizations do become important.  
            Above mentioned equation is used for determining the SOFC voltage (in fact for fuel cell voltage in general). 
            This approach results in good agreement with particular experimental data (for which adequate factors were 
            obtained) and poor agreement for other than original experimental working parameters. Moreover, most of the 
            equations used require the addition of numerous factors which are difficult or impossible to determine. 
            It makes very difficult any optimizing process of the SOFC working parameters as well as design 
            architecture configuration selection. Because of those circumstances a few other equations were proposed:
            '''
            st.latex(
                r'''E_{SOFC}=\frac{E_{max}-i_{max} \cdot \eta_f\cdot r_1}{\frac{r_1}{r_2}\cdot (1-\eta_f)+1}''')
            '''
            where:  
            >$E_{SOFC}=$ cell voltage  
            >$E_{max}=$ maximum voltage given by the Nernst equation  
            >$i_{max}=$ maximum current density (for given fuel flow)  
            >$\eta_f=$ fuel utilization factor  
            >$r_1=$ ionic specific resistance of the electrolyte  
            >$r_2=$ electric specific resistance of the electrolyte.  
            '''
            '''
            This method was validated and found to be suitable for optimization and sensitivity studies in plant-level 
            modelling of various systems with solid oxide fuel cells. With this mathematical description it is 
            possible to account for different properties of the SOFC. There are many parameters which impact cell working 
            conditions, e.g. electrolyte material, electrolyte thickness, cell temperature, inlet and outlet gas compositions 
            at anode and cathode, and electrode porosity, just to name some. The flow in these systems is often calculated 
            using the [Navier–Stokes equations](https://en.wikipedia.org/wiki/Navier–Stokes_equations).
            '''
            with st.expander('Ohmic polarization'):
                '''
                Ohmic losses in an SOFC result from ionic conductivity through the electrolyte and electrical resistance 
                offered to the flow of electrons in the external electrical circuit. This is inherently a materials property 
                of the crystal structure and atoms involved. However, to maximize the ionic conductivity, several methods 
                can be done. Firstly, operating at higher temperatures can significantly decrease these ohmic losses. 
                Substitutional doping methods to further refine the crystal structure and control defect concentrations 
                can also play a significant role in increasing the conductivity. Another way to decrease ohmic resistance 
                is to decrease the thickness of the electrolyte layer.
                '''
                '''
                An ionic specific resistance of the electrolyte as a function of temperature can be described by the following relationship:
                '''
                st.latex(r''' r_1=\frac{\delta}{\sigma}''')
                '''
                where:  
                >$\delta$ - electrolyte thickness  
                >$\sigma$ - ionic conductivity  
                '''
                '''
                The ionic conductivity of the solid oxide is defined as follows:
                '''
                st.latex(r''' \sigma = \sigma_0 \cdot e^{\frac{-E}{R \cdot T}}''')
                '''
                where:  
                >$\sigma_{0}$ and $E$ – factors depended on electrolyte materials  
                >$T$ – electrolyte temperature  
                >$R$ – ideal gas constant  
                '''
            with st.expander('Concentration polarization'):
                '''
                The concentration polarization is the result of practical limitations on mass transport within the cell and 
                represents the voltage loss due to spatial variations in reactant concentration at the chemically active sites. 
                This situation can be caused when the reactants are consumed by the electrochemical reaction faster than they 
                can diffuse into the porous electrode, and can also be caused by variation in bulk flow composition. The latter 
                is due to the fact that the consumption of reacting species in the reactant flows causes a drop in reactant 
                oncentration as it travels along the cell, which causes a drop in the local potential near the tail end of the cell.
                '''
                '''
                The concentration polarization occurs in both the anode and cathode. The anode can be particularly problematic, 
                as the oxidation of the hydrogen produces steam, which further dilutes the fuel stream as it travels along the 
                length of the cell. This polarization can be mitigated by reducing the reactant utilization fraction or increasing 
                the electrode porosity, but these approaches each have significant design trade-offs.
                '''
            with st.expander('Activation polarization'):
                '''
                The activation polarization is the result of the kinetics involved with the electrochemical reactions. 
                Each reaction has a certain activation barrier that must be overcome in order to proceed and this barrier 
                leads to the polarization. The activation barrier is the result of many complex electrochemical reaction steps 
                where typically the rate limiting step is responsible for the polarization. The polarization equation shown below 
                is found by solving the [Butler–Volmer](https://en.wikipedia.org/wiki/Butler–Volmer_equation) equation in the 
                high current density regime (where the cell typically operates), and can be used to estimate the activation polarization:
                '''
                st.latex(
                    r''' \eta_{act}=\frac{RT}{\beta zF} \cdot ln \left( \frac{i}{i_{0}} \right)''')
                '''
                where:  
                >$R =$ gas constant  
                >${T}_0 =$ operating temperature  
                >$z =$ electrons associated with the electrochemical reaction  
                >$F =$ Faraday's constant  
                >$i =$ operating current  
                >$i_{0} =$ exchange current density  
                '''
                '''
                The polarization can be modified by microstructural optimization. The Triple Phase Boundary (TPB) length, 
                which is the length where porous, ionic and electronically conducting pathways all meet, directly relates to 
                the electrochemically active length in the cell. The larger the length, the more reactions can occur and thus 
                the less the activation polarization. Optimization of TPB length can be done by processing conditions to affect 
                microstructure or by materials selection to use a mixed ionic/electronic conductor to further increase TPB length.
                '''
            '''
            ***
            '''

        with whatnext_container:
            '''
            ### What next?
            Based on the above equations and some more machine learning algorithms we 
            have prepared a mathematical model that describes SOFC single cell. Please, 
            follow up all the works for better understanding of SOFC theoretical basics.
            '''
            '''
            To proceed, select the next section in the left sidebar. 
            '''

        with st.container():
            '''
            ### Introduction & theoretical background  
            A **current–voltage characteristic** or **I–V curve** 
            (current–voltage curve) is a relationship, typically represented 
            as a chart or graph, between the electric current through a circuit, 
            device, or material, and the corresponding voltage, or potential 
            difference across it.
            '''
            '''
            A single SOFC cell operated with hydrogen and oxygen provides at 
            equilibrium a theoretical reversible (Nernst) or open circuit voltage (OCV) 
            of $1.229 V$ at standard conditions $(T = 273.15 K, p = 1 atm$. 
            With the standard electrode potential $E_0$, universal gas constant $R$, 
            temperature $T$, Faraday’s constant $F$, molar concentration $x$ and pressure $p$, 
            the OCV is given by the following:
            '''
            st.latex(
                r'''E(p,T)=E^0(p^0,T)-\frac{RT}{2F}\ln \left( \frac{x_{H_2O}}{x_{H_2} \cdot \sqrt{x_{O_2}}} \right) + \frac{RT}{4F}\ln{p}''')
            '''
            However, the actual measured OCV will often fall slightly below $E$ by $U_L$ which 
            represents losses in potential due to residual electronic conduction in the electrolyte 
            and possibly also cross over of gases via micro cracks and fissures in the electrolyte. 
            The thermodynamically obtainable cell voltage Vcell at OCV, moreover, depends on the 
            used fuel and particularly, on operation temperature and pressure. For example, the 
            OCV of an atmospheric SOFC operating on hydrogen and oxygen is about $0.908 V$ at 
            $1000 \N{DEGREE SIGN}C$. The actual voltage is lower than the theoretical model due 
            to reaction, charge and mass transfer losses. For more information, see the 
            **Introduction** chapter, **Polarizations** section. As shown in figure below, 
            the performance of SOFC cell can be illustrated using a polarization curve that can 
            be broken into three sections:  
            (1) activation losses,  
            (2) ohmic losses,   
            (3) mass transport losses.  
            Therefore, the operating voltage of the cell can be represented as the departure 
            from ideal voltage caused by these three losses (polarizations).
            '''
            st.markdown('<p style="text-align: center; "><img src="https://raw.githubusercontent.com/MelnikovAP/sofc-streamlit/master/Images/polarization_curve.png" width=500/></p>',
                        unsafe_allow_html=True)

if tab_selected == works[0] and lang_selected == langs[0]:
    with st.container():
        st.markdown('''
        <div style="text-align: justify;">
        Современная энергетика базируется в основном на электромеханическом способе преобразования энергии, когда тепловая энергия,
         выделяющаяся при сгорании топлива сначала преобразуется в механическую (обычно вращения), которая в свою очередь 
         преобразуется в электрическую посредством электрогенератора. Однако, обусловленный развитием цивилизации рост 
         энергопотребления при исчерпаемости ископаемых энергоносителей, возрастание их стоимости и близкая к предельной 
         экологическая нагрузка побуждают человечество предпринимать усилия по повышению эффективности преобразования энергии 
         первичных источников в электрическую и развивать альтернативные способы ее производства. Такой альтернативой могут
          стать работающие по принципу прямого преобразования энергии топливные элементы (ТЭ), которые позволяют сразу получать из 
          энергии химических связей топлива электрическую без промежуточного перехода в механическую энергию. Такой процесс 
          получения электроэнергии в ТЭ значительно более эффективен, чем в традиционно используемых в энергетике электромеханических 
          преобразователях. В ТЭ нет движущихся частей что значительно увеличивает КПД процесса.
        <p></p></div>
        ''', unsafe_allow_html=True)

        st.markdown('''
        <div style="text-align: justify;">
        Нельзя сказать, что ТЭ уже являются обыденными источниками энергии, но несомненно, технологии ТЭ переживают бурное развитие. 
        ТЭ уже находят применение в стационарных энергоустановках широкого диапазона мощностей, транспортных силовых установках, 
        переносных и портативных источниках электропитания.
        <p></p></div>
        ''', unsafe_allow_html=True)

        st.markdown('''
        <div style="text-align: justify;">
        В ТЭ химическая энергия топлива непосредственно преобразуется в электричество в процессе бесшумной и беспламенной 
        электрохимической реакции. Электролит различных ТЭ может находиться в твёрдом (полимеры, гибридные материалы, керамика)
         или жидком (раствор или расплав) агрегатном состоянии и должен обладать высокой ионной проводимостью (<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline-block" title="O^{2-} "> <mrow> <msup> <mrow> O </mrow> <mrow> <mn>2</mn> <mo>-</mo> </mrow> </msup> </mrow></math>, <math xmlns="http://www.w3.org/1998/Math/MathML" display="inline-block" title="H^{+} "> <mrow> <msup> <mrow> H </mrow> <mrow> <mo>+</mo> </mrow> </msup> </mrow></math> 
         и другие ионы, в зависимости от типа ТЭ) в сочетании с пренебрежимо малой электронной проводимостью. Тип электролита 
         часто служит основой при классификации ТЭ.
        <p></p></div>
        ''', unsafe_allow_html=True)

        st.markdown('''
        <div style="text-align: justify;">
        Среди различных типов топливных элементов, для энергоустановок стационарного назначения и целого ряда транспортных приложений 
        наиболее подходящими представляются высокотемпературные твердооксидные топливные элементы (ТОТЭ).
        <p></p></div>
        ''', unsafe_allow_html=True)

        st.markdown('''
        <div style="text-align: justify;">
        ТОТЭ являются топливными элементами с самой высокой рабочей температурой. Рабочая температура ТОТЭ в общем случае может 
        варьироваться от 600°C до 900°C. КПД производства электрической энергии у ТОТЭ – самый высокий из всех топливных элементов и в 
        принципе может достигать 70%.
        <p></p></div>
        ''', unsafe_allow_html=True)


        st.markdown('''
        <div style="text-align: justify;">
        Важнейшими преимуществами ТОТЭ являются широкий спектр потребляемых видов топлива – газообразные и жидкие углеводороды, спирты.
        <p></p></div>
        ''', unsafe_allow_html=True)


        st.markdown('''
        <div style="text-align: justify;">
        Основными компонентами любого единичного ТЭ являются катод, анод и разделяющий их электролит.
        <p></p></div>
        ''', unsafe_allow_html=True)

        st.markdown('''
        <div style="text-align: justify;">
        Принципиальная рабочая схема любого ТЭ выглядит довольно просто. К одной стороне ТЭ – аноду – необходимо подавать топливо (синтез-газ, 
        содержащий водород и окись углерода, для ТОТЭ; чистый водород для наиболее распространённых низкотемпературных видов ТЭ или 
        другие виды топлив для отдельных типов ТЭ), к другой стороне – катоду – подают воздух (или чистый кислород).
        <p></p></div>
        ''', unsafe_allow_html=True)

        st.markdown('''
        <div style="text-align: justify;">
        В ТОТЭ анод и катод представляют собой тонкие слои керамических и металлокерамических композитных материалов разного 
        состава с открытой пористостью. Пористая структура электродов необходима, т.к. именно на развитой поверхности 
        многочисленных пор происходят основные химические реакции.
        <p></p></div>
        ''', unsafe_allow_html=True)    

        st.markdown('''
        <div style="text-align: justify;">
        Отличительной особенностью ТОТЭ является то, что для работы при высоких температурах используемый электролит представляет 
        собой тонкую твердую керамическую структуру на основе оксидов металлов, часто в составе содержащих 
        иттрий и цирконий, которая является проводником ионов кислорода (<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline-block" title="O^{2-} "> <mrow> <msup> <mrow> O </mrow> <mrow> <mn>2</mn> <mo>-</mo> </mrow> </msup> </mrow></math>). Именно благодаря особенностям электролита 
        твёрдооксидный топливный элемент получил своё название.  
        <p></p></div>
        ''', unsafe_allow_html=True)    

        st.markdown('''
        <div style="text-align: justify;">
        На рисунке 1 показана упрощённая принципиальная схема ТОТЭ.  
        <p></p></div>
        ''', unsafe_allow_html=True) 

        col1, col2, col3 = st.columns([1,6,1])
        col2.image(load_image(os.path.dirname(__file__) + r"\Images\intro_fig1.jpg"), width=500, caption="Рисунок 1. − Принципиальная схема ТОТЭ")

        st.markdown('''
        <div style="text-align: justify;">
        Рассмотрим электрохимические реакции в ТОТЭ на простом примере водородного топлива. На катоде протекает реакция восстановления кислорода:
        <math xmlns="http://www.w3.org/1998/Math/MathML" display="block" title="O_2 + 4\overline{e} = 2O^{2-} "> <mrow> <msub> <mrow> O </mrow> <mrow> <mn>2</mn> </mrow> </msub> <mo>+</mo> <mn>4</mn> <mover accent="true"> <mrow> <mi>e</mi> </mrow> <mo>¯</mo> </mover> <mo>=</mo> <mn>2</mn> <msup> <mrow> O </mrow> <mrow> <mn>2</mn> <mo>-</mo> </mrow> </msup> </mrow>,</math>  
        <p></p></div>
        ''', unsafe_allow_html=True)    

        st.markdown('''
        <div style="text-align: justify;">
        а на аноде – окисления топлива:
        <math xmlns="http://www.w3.org/1998/Math/MathML" display="block" title="H_2 + O^{2-} = H_2O+ 2\overline{e} "> <mrow> <msub> <mrow> H </mrow> <mrow> <mn>2</mn> </mrow> </msub> <mo>+</mo> <msup> <mrow> O </mrow> <mrow> <mn>2</mn> <mo>-</mo> </mrow> </msup> <mo>=</mo> <msub> <mrow> H </mrow> <mrow> <mn>2</mn> </mrow> </msub> O <mo>+</mo> <mn>2</mn> <mover accent="true"> <mrow> <mi>e</mi> </mrow> <mo>¯</mo> </mover> </mrow></math>  
        <p></p></div>
        ''', unsafe_allow_html=True)   

        st.markdown('''
        <div style="text-align: justify;">
        Ионы кислорода движутся через ионопроводящий твёрдый электролит от катода к аноду, где соединяются с водородом. 
        Продуктом реакции в этом случае является вода. Общая реакция окисления:
        <math xmlns="http://www.w3.org/1998/Math/MathML" display="block" title="2H_2 + O_2 \rightarrow 2H_2O "> <mrow> <mn>2</mn> <msub> <mrow> H </mrow> <mrow> <mn>2</mn> </mrow> </msub> <mo>+</mo> <msub> <mrow> O </mrow> <mrow> <mn>2</mn> </mrow> </msub> <mo>→</mo> <mn>2</mn> <msub> <mrow> H </mrow> <mrow> <mn>2</mn> </mrow> </msub> O </mrow>,</math>
        <p></p></div>
        ''', unsafe_allow_html=True)   

        st.markdown('''
        <div style="text-align: justify;">
        такая же, как и при горении водорода. Однако в ТЭ потоки топлива и окислителя не смешиваются, а реакции окисления топлива 
        и восстановления кислорода, как и в батарейках, пространственно разделены и проходят на разных электродах – соответственно, 
        процесс «сжигания» протекает, только если элемент попутно выдает ток во внешнюю цепь, вырабатывая электричество.  
        <p></p></div>
        ''', unsafe_allow_html=True)   

        st.markdown('''
        <div style="text-align: justify;">
        Электролит ТОТЭ обеспечивает транспорт ионов кислорода <math xmlns="http://www.w3.org/1998/Math/MathML" display="inline-block" title="O^{2-} "> <mrow> <msup> <mrow> O </mrow> <mrow> <mn>2</mn> <mo>-</mo> </mrow> </msup> </mrow></math> от катода к аноду и разделяет два газовых объёма: 
        топливный и окислительный. Высоко избирательные ионопроводящие свойства твёрдого электролита (при минимально 
        возможной электронной проводимости) достигаются тщательным подбором состава (соотношения основных компонентов 
        и легирующих добавок) и специальными методами изготовления из него прецизионных тонкоплёночных структур, 
        обеспечивающих газовую непроницаемость и высокую ионную проводимость. Ионопроводящая эффективность электролита 
        напрямую зависит от толщины и температуры – падает с ростом первой и растёт с ростом второй. Чем тоньше можно 
        сделать плёнку электролита, тем лучшую ионную проводимость он обеспечивает. При этом очень высокая температура 
        процесса в ТОТЭ помимо преимущества в виде высокой эффективности имеет и существенные недостатки, а именно, 
        необходимость применения в ТОТЭ дорогостоящих материалов, способных выдерживать такие критические условия на 
        протяжении длительного времени. Уменьшая толщину электролита, можно в некоторой степени снизить температуру 
        процесса в ТОТЭ, что очень важно с точки зрения подбора более доступных и дешёвых материалов для ТОТЭ.  
        <p></p></div>
        ''', unsafe_allow_html=True)   

        st.markdown('''
        <div style="text-align: justify;">
        Не менее важным является подбор эффективных материалов для катода и анода, т.к. они должны отвечать целому 
        комплексу специфических, часто противоречивых требований: химическая стойкость и стабильность при высоких 
        температурах в условиях характерных рабочих сред, высокая проницаемость для прохождения рабочих газовых сред 
        при достаточной прочности, хорошая адгезия и взаимная совместимость по коэффициенту термического расширения с 
        материалом электролита в широком диапазоне температур, определённые каталитические свойства анодного материала, 
        интенсифицирующие проведение целевой электрохимической реакции.
        <p></p></div>
        ''', unsafe_allow_html=True)   

        st.markdown('''
        <div style="text-align: justify;">
        Следует отметить для представления сложностей в производстве, что толщина электролита, в зависимости от 
        подхода производителей, варьируется примерно от 10 до 150 микрон, то есть от 0.01 до 0.15 мм. Толщины анода 
        и катода имеют схожий порядок величин.
        <p></p></div>
        ''', unsafe_allow_html=True)  

        st.markdown('''
        <div style="text-align: justify;">
        Батареи топливных элементов, которые находят непосредственное применение в электрохимических генераторах, 
        состоят из определённого количества собранных совместно единичных ТЭ (трубок) и других вспомогательных элементов, 
        позволяющих работать отдельным ТЭ совместно, обеспечивающих равномерную подачу топлива и окислителя к анодам и 
        катодам всех ТЭ одновременно, вывод отработанных продуктов реакции и электрическую коммутацию, внутреннюю и внешнюю.
        <p></p></div>
        ''', unsafe_allow_html=True)  

        st.markdown('''
        <div style="text-align: justify;">
        Рисунок 2 представляет примерное устройство условного единичного блока планарного ТОТЭ, с организацией подвода 
        топлива и окислителя к рабочим поверхностям электродов.
        <p></p></div>
        ''', unsafe_allow_html=True)  

        col1, col2, col3 = st.columns([1,6,1])
        col2.image(load_image(os.path.dirname(__file__) + r"\Images\intro_fig2.jpg"), width=500, caption="Рисунок 2. − Пример устройства единичного планарного ТОТЭ")


        st.markdown('''
        <div style="text-align: justify;">
        Если изображённую на рисунке 2 единичную структуру многократно повторить при последовательном 
        её сложении по оси У, получится батарея планарных топливных элементов, как показано на рисунке 3.
        <p></p></div>
        ''', unsafe_allow_html=True)  

        col1, col2, col3 = st.columns([1,6,1])
        col2.image(load_image(os.path.dirname(__file__) + r"\Images\intro_fig3.jpg"), width=500, caption="Рисунок 3. − Пример формирования планарной батареи ТОТЭ")


        st.markdown('''
        <div style="text-align: justify;">
        Батарея трубчатых ТОТЭ представляет собой общий коллектор из множества параллельных трубок, внутренняя полость которых предназначена для движения одной 
        из рабочих сред (топлива или окислителя), а наружное пространство – для другой, как в упрощённом виде показано на рисунке 4.
        <p></p></div>
        ''', unsafe_allow_html=True)  

        col1, col2, col3 = st.columns([1,6,1])
        col2.image(load_image(os.path.dirname(__file__) + r"\Images\intro_fig4.jpg"), width=500, caption="Рисунок 4. − Пример организации батареи трубчатых ТОТЭ")


# @st.cache(allow_output_mutation=True)
class VI():
    def makePlotVI(self,params=None):
        temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc = params
        _result = get_ivc(temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc)
        _p = figure(
            x_axis_label=td[lang_selected][18] + ', ' + td[lang_selected][9],
            y_axis_label=td[lang_selected][17] + ', '+ td[lang_selected][15],
            plot_width=500, plot_height=300,
        )
        _p.sizing_mode = "scale_both"
        _p.line(_result['i'], (_result['Eload']-_result['eta_con']),  line_width=2)
        _p.axis.axis_label_text_font_style = 'bold'
        # _p.add_layout(_text)
        return _p


class PI():
    def makePlotPI(self,params=None):
        temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc = params
        _result = get_ivc(temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc)
        _p = figure(
            x_axis_label=td[lang_selected][18] + ', ' + td[lang_selected][9],
            y_axis_label=td[lang_selected][19] + ', '+ td[lang_selected][16],
            plot_width=500, plot_height=300,
        )
        _p.sizing_mode = "scale_both"
        _p.line(_result['i'], _result['i'] * (_result['Eload']-_result['eta_con']),  line_width=2)
        _p.axis.axis_label_text_font_style = 'bold'
        # _p.add_layout(_text)
        return _p


for k in ["vi1","vi2","vi3","vi4","vi5","vi6","vi7","vi8","pi1","pi2"]:
    if k not in st.session_state.keys():
        st.session_state[k] = [
            temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
# if "tvi2" not in st.session_state.keys():
#     st.session_state['tvi2'] = [
#         temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
# if "tvi3" not in st.session_state.keys():
#     st.session_state['tvi3'] = [
#         temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
# if "tvi4" not in st.session_state.keys():
#     st.session_state['tvi4'] = [
#         temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
# if "svi1" not in st.session_state.keys():
#     st.session_state['svi1'] = [
#         temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
# if "svi2" not in st.session_state.keys():
#     st.session_state['svi2'] = [
#         temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
# if "svi3" not in st.session_state.keys():
#     st.session_state['svi3'] = [
#         temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
# if "svi4" not in st.session_state.keys():
#     st.session_state['svi4'] = [
#         temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
# if "mvi1" not in st.session_state.keys():
#     st.session_state['mvi1'] = [
#         temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
# if "mvi2" not in st.session_state.keys():
#     st.session_state['mvi2'] = [
#         temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]



task_title = lambda tsk,ttl: '''
    <div style="text-align: center; padding:70px 0px 10px 0px"> 
    <span style="padding:2px 7px 2px 7px;color:white;background-color: rgb(77,77,77);border-radius:3px;display:inline-flex;justify-content:center;align-items:center;font-weight:bold;">{0:s}</span><span style="padding-left:10px;">{1:s}</span>
    <p></p></div>
    '''.format(tsk,ttl)

task_subtitle = lambda ttl: '''
    <div style="text-align: center;font-weight:bold;"> 
    {0:s}
    <p></p></div>
    '''.format(ttl)

plt_statetitle = lambda sts: '''
        <div style="text-align: right;font-size: 12px; padding-top:2px;">        
        <b>Saved state:</b> {0:.0f} K, {1:.1f} atm, {2:.2f}⋅10⁻⁶ S/m, {3:.2f} µm, {4:.2f} mA/cm²
        <p></p></div>
        '''.format(*st.session_state[sts])



plt_vi1 = VI()
plt_vi2 = VI()
plt_vi3 = VI()
plt_vi4 = VI()
plt_vi5 = VI()
plt_vi6 = VI()
plt_vi7 = VI()
plt_vi8 = VI()
plt_pi1 = PI()
plt_pi2 = PI()

if tab_selected == works[1] and lang_selected == langs[0]:

    with st.container():
        st.header("Effect of temperature on the performance of fuel cell")
        st.markdown('''
        <div style="text-align: justify;">        
        The ideal performance of a fuel cell depends on the electrochemical reactions that occur with different fuels and oxygen.
        Standard-state reversible fuel cell voltages (E₀ values) are only useful under standard-state conditions (room temperature, atmospheric pressure, unit activities of all species). 
        Fuel cells are frequently operated under conditions that vary greatly from the standard state. For example, high-temperature fuel cells operate at 700–1000°C, 
        automotive fuel cells often operate under 3–5 atm of pressure, and almost all fuel cells cope with variations in the concentration (and therefore activity) of reactant species. 
        In the following sections, we systematically define how reversible fuel cell voltages are affected by departures from the standard state. 
        First, the influence of temperature on the reversible fuel cell voltage will be explored, then the influence of pressure.
        <p></p></div>
        ''', unsafe_allow_html=True)


       
        with st.container():
            
            st.markdown('''
            To understand how the reversible voltage varies with temperature, lets look at the differential expression for the Gibbs free energy $G$ - one of
            the most important themodynamic relation for chemical reactions in the fuel cell:
            ''')

            st.latex(r'''
            \tag{1} dG = −SdT + Vdp
            ''')

            st.markdown('''
            where $S$ is entropy, $T$, $V$ and $p$ are thermodynamic parameters.  

            From this equation, we find the relation between molar reaction quantities $\Delta g$  and $\Delta s$:
            ''')

            st.latex('''
            (d (\Delta g)/dT)_p = −\Delta s
            ''')

            st.markdown('''
            The Gibbs free energy $\Delta g$ is related to the cell voltage $E$, Faraday constant $F$  
            and the number of moles of transferred electrons $n$ by the equation:
            ''')

            st.latex(r'''
            \tag{2} \Delta g = -nFE
            ''')   

            st.markdown('''
            When we combine all previous expressions we can find how a fuel cell voltage $E$ varies as a function of temperature:
            ''')         

            st.latex(r'''
            (dE/dT)_p = \frac{\Delta s}{nF}
            ''')   


            st.markdown('''
            or in intagrated form:
            ''')

            st.latex(r'''
            \tag{3} E_T=E_0+ nF \cdot \Delta s \cdot (T−T_0) 
            ''')   

            st.markdown('''
            where  $E_T$ is the reversible cell voltage at an arbitrary temperature $T$ and constant pressure. 
            Generally, we assume $\Delta s$ to be independent of temperature. 
            If a more accurate value of $E_T$ is required, it may be calculated by integrating the heat-capacity-related temperaturedependence of $\Delta s$.
            ''', unsafe_allow_html=True)

            st.markdown('''
            A task below wil help you to understand the influence of working temperature on the fuel cell voltage.
            ''')


            st.markdown(task_title("Задание 1.1","Области низких и высоких температур"), unsafe_allow_html=True)

            st.markdown(task_subtitle("Работа с моделью"), unsafe_allow_html=True)

            st.markdown('''
            * Using the sidebars at the `Model parameters` panel change cell temperature to 900 K. 
            * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
            * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.
            ''')




        
            col1, col2 = st.columns([2,1])
            chk_vi1 = col2.checkbox(td[lang_selected][14], key="chk_vi1")
            if chk_vi1 and btn_runsimulation:
                st.session_state["vi1"] = [
                    temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
            st.bokeh_chart(plt_vi1.makePlotVI(st.session_state["vi1"]))
            col1.markdown(plt_statetitle("vi1"), unsafe_allow_html=True)


            st.markdown('''
            * Using the sidebars at the `Model parameters` panel change cell temperature to 1000 K. 
            * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
            * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.
            ''')


            col1, col2 = st.columns([2,1])
            chk_vi2 = col2.checkbox(td[lang_selected][14], key="chk_vi2")
            if chk_vi2 and btn_runsimulation:
                st.session_state["vi2"] = [
                    temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]

            st.bokeh_chart(plt_vi2.makePlotVI(st.session_state["vi2"]))
            col1.markdown(plt_statetitle("vi2"), unsafe_allow_html=True)

            st.markdown(task_subtitle("Анализ результатов"), unsafe_allow_html=True)


            st.markdown('''
            Compare the obtained Voltage-Current curves  and answer the questions:
             
            * How the fuel cell voltage $V$ reacts on temperature changing?
            * How your assumptions relate with the above linear expression (3)? 

            ''')
        with st.container():
            st.markdown('''
                Like temperature effects, the pressure effects on cell voltage may also be calculated starting from the
                differential expression for the Gibbs free energy (1) and its electrochemical definition (2): 
                ''')

            st.latex(r'''
                (d G/dp)_T = \Delta V \Rightarrow \left(\frac{d (\Delta g)}{dp}\right)_T = \Delta v 
                ''')

            st.markdown('''
                and, finally:
                ''')

            st.latex(r'''
                \tag{4} (dE/dp)_T = -\frac{\Delta v}{nF} 
                ''')

            st.markdown('''
                Expression (4) means the variation of the reversible cell voltage with pressure is related to the volume change of the reaction. 
                If the volume change of the reaction is negative (if fewer moles of gas are generated by the reaction than consumed, for instance), 
                then the cell voltage will increase with increasing pressure. This is an example of *Le Chatelier’s principle*: 
                Increasing the pressure of the system favors the reaction direction that relieves the stress on the system.
                ''')

            st.markdown('''
                To understand how the fuel cell  works at non-ambient pressure, complete the following tasks:
                ''')


            st.markdown(task_title("Задание 1.2","Области нормального и высокого давлений"), unsafe_allow_html=True)

            st.markdown(task_subtitle("Работа с моделью"), unsafe_allow_html=True)

            st.markdown('''
            * Using the sidebars at the `Model parameters` panel change acting pressure value to 1.1 atm.   
            * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
            * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.
            ''')

        
            col1, col2 = st.columns([2,1])
            chk_vi3 = col2.checkbox(td[lang_selected][14], key="chk_vi3")
            if chk_vi3 and btn_runsimulation:
                st.session_state["vi3"] = [
                    temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
            st.bokeh_chart(plt_vi3.makePlotVI(st.session_state["vi3"]))
            col1.markdown(plt_statetitle("vi3"), unsafe_allow_html=True)
    

            st.markdown('''
            * Using the sidebars at the `Model parameters` panel change acting pressure value to 1.8 atm. 
            * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
            * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.
            ''')

        
            col1, col2 = st.columns([2,1])
            chk_vi4 = col2.checkbox(td[lang_selected][14], key="chk_vi4")
            if chk_vi4 and btn_runsimulation:
                st.session_state["vi4"] = [
                    temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]

            st.bokeh_chart(plt_vi4.makePlotVI(st.session_state["vi4"]))
            col1.markdown(plt_statetitle("vi4"), unsafe_allow_html=True)





            st.markdown(task_subtitle("Анализ результатов"), unsafe_allow_html=True)

            st.markdown('''
            Compare the obtained Voltage-Current curves and answer the questions:
            
            * How the fuel cell voltage $V$ reacts on the pressure changing?
            * How do you think, is $V(p)$  follows *Le Chatelier’s principle*? 

            ''')


if tab_selected == works[1] and lang_selected == langs[1]:

    cflag = False

    if cflag:
        st.markdown('''
        Under construction...
        ''')
    else:

        with st.container():
            st.header(
                "Effect of temperature on the performance of fuel cell")
            st.markdown('''
            <div style="text-align: justify;">
            The ideal performance of a fuel cell depends on the electrochemical reactions that occur with different fuels and oxygen.
            Standard-state reversible fuel cell voltages (E₀ values) are only useful under standard-state conditions (room temperature, atmospheric pressure, unit activities of all species).
            Fuel cells are frequently operated under conditions that vary greatly from the standard state. For example, high-temperature fuel cells operate at 700–1000°C,
            automotive fuel cells often operate under 3–5 atm of pressure, and almost all fuel cells cope with variations in the concentration (and therefore activity) of reactant species.
            In the following sections, we systematically define how reversible fuel cell voltages are affected by departures from the standard state.
            First, the influence of temperature on the reversible fuel cell voltage will be explored, then the influence of pressure.
            <p></p></div>
            ''', unsafe_allow_html=True)



            with st.container():
                
                st.markdown('''
                To understand how the reversible voltage varies with temperature, lets look at the differential expression for the Gibbs free energy $G$ - one of
                the most important themodynamic relation for chemical reactions in the fuel cell:
                ''')

                st.latex(r'''
                \tag{1} dG = −SdT + Vdp
                ''')

                st.markdown('''
                where $S$ is entropy, $T$, $V$ and $p$ are thermodynamic parameters.  

                From this equation, we find the relation between molar reaction quantities $\Delta g$  and $\Delta s$:
                ''')

                st.latex('''
                (d (\Delta g)/dT)_p = −\Delta s
                ''')

                st.markdown('''
                The Gibbs free energy $\Delta g$ is related to the cell voltage $E$, Faraday constant $F$  
                and the number of moles of transferred electrons $n$ by the equation:
                ''')

                st.latex(r'''
                \tag{2} \Delta g = -nFE
                ''')   

                st.markdown('''
                When we combine all previous expressions we can find how a fuel cell voltage $E$ varies as a function of temperature:
                ''')         

                st.latex(r'''
                (dE/dT)_p = \frac{\Delta s}{nF}
                ''')   


                st.markdown('''
                or in intagrated form:
                ''')

                st.latex(r'''
                \tag{3} E_T=E_0+ nF \cdot \Delta s \cdot (T−T_0) 
                ''')   

                st.markdown('''
                where  $E_T$ is the reversible cell voltage at an arbitrary temperature $T$ and constant pressure. 
                Generally, we assume $\Delta s$ to be independent of temperature. 
                If a more accurate value of $E_T$ is required, it may be calculated by integrating the heat-capacity-related temperaturedependence of $\Delta s$.
                ''', unsafe_allow_html=True)

                st.markdown('''
                A task below wil help you to understand the influence of working temperature on the fuel cell voltage.
                ''')


                st.markdown(task_title("Task 1.1","Low and high temperatures"), unsafe_allow_html=True)

                st.markdown(task_subtitle("Simulations"), unsafe_allow_html=True)

                st.markdown('''
                * Using the sidebars at the `Model parameters` panel change cell temperature to 900 K. 
                * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
                * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.
                ''')




            
                col1, col2 = st.columns([2,1])
                chk_vi1 = col2.checkbox(td[lang_selected][14], key="chk_vi1")
                if chk_vi1 and btn_runsimulation:
                    st.session_state["vi1"] = [
                        temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
                st.bokeh_chart(plt_vi1.makePlotVI(st.session_state["vi1"]))
                col1.markdown(plt_statetitle("vi1"), unsafe_allow_html=True)


                st.markdown('''
                * Using the sidebars at the `Model parameters` panel change cell temperature to 1100 K. 
                * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
                * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.
                ''')


                col1, col2 = st.columns([2,1])
                chk_vi2 = col2.checkbox(td[lang_selected][14], key="chk_vi2")
                if chk_vi2 and btn_runsimulation:
                    st.session_state["vi2"] = [
                        temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]

                st.bokeh_chart(plt_vi2.makePlotVI(st.session_state["vi2"]))
                col1.markdown(plt_statetitle("vi2"), unsafe_allow_html=True)

                st.markdown(task_subtitle("Results analisys"), unsafe_allow_html=True)


                st.markdown('''
                Compare the obtained Voltage-Current curves  and answer the questions:
                
                * How the fuel cell voltage $V$ reacts on temperature changing?
                * Does it follow from the above linear expression (3)? 

                ''')
            with st.container():
                st.markdown('''
                    Like temperature effects, the pressure effects on cell voltage may also be calculated starting from the
                    differential expression for the Gibbs free energy (1) and its electrochemical definition (2): 
                    ''')

                st.latex(r'''
                    (d G/dp)_T = \Delta V \Rightarrow \left(\frac{d (\Delta g)}{dp}\right)_T = \Delta v 
                    ''')

                st.markdown('''
                    and, finally:
                    ''')

                st.latex(r'''
                    \tag{4} (dE/dp)_T = -\frac{\Delta v}{nF} 
                    ''')

                st.markdown('''
                    Expression (4) means the variation of the reversible cell voltage with pressure is related to the volume change of the reaction. 
                    If the volume change of the reaction is negative (if fewer moles of gas are generated by the reaction than consumed, for instance), 
                    then the cell voltage will increase with increasing pressure. This is an example of *Le Chatelier’s principle*: 
                    Increasing the pressure of the system favors the reaction direction that relieves the stress on the system.
                    ''')

                st.markdown('''
                    To understand how the fuel cell  works at non-ambient pressure, complete the following tasks:
                    ''')


                st.markdown(task_title("Task 1.2","Normal and high pressures"), unsafe_allow_html=True)

                st.markdown(task_subtitle("Simulations"), unsafe_allow_html=True)

                st.markdown('''
                * Using the sidebars at the `Model parameters` panel change acting pressure value to 1.1 atm. 
                * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
                * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.
                ''')

            
                col1, col2 = st.columns([2,1])
                chk_vi3 = col2.checkbox(td[lang_selected][14], key="chk_vi3")
                if chk_vi3 and btn_runsimulation:
                    st.session_state["vi3"] = [
                        temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
                st.bokeh_chart(plt_vi3.makePlotVI(st.session_state["vi3"]))
                col1.markdown(plt_statetitle("vi3"), unsafe_allow_html=True)
        

                st.markdown('''
                * Using the sidebars at the `Model parameters` panel change acting pressure value to 1.8 atm.  
                * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
                * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.
                ''')

            
                col1, col2 = st.columns([2,1])
                chk_vi4 = col2.checkbox(td[lang_selected][14], key="chk_vi4")
                if chk_vi4 and btn_runsimulation:
                    st.session_state["vi4"] = [
                        temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]

                st.bokeh_chart(plt_vi4.makePlotVI(st.session_state["vi4"]))
                col1.markdown(plt_statetitle("vi4"), unsafe_allow_html=True)





                st.markdown(task_subtitle("Анализ результатов"), unsafe_allow_html=True)

                st.markdown('''
                Compare the obtained Voltage-Current curves and answer the questions:
                
                * How the fuel cell voltage $V$ reacts on the pressure changing?
                * Does it follow from the above linear expression (4)? 

                ''')


if tab_selected == works[2] and lang_selected == langs[0]:

    with st.container():

        st.header("Ohmic losses")

        st.markdown(r'''
        Ohmic losses ($E_{\text{Ohm}}$) in an SOFC result from ionic conductivity through the electrolyte and electrical resistance 
        offered to the flow of electrons in the external electrical circuit. This is inherently a materials property 
        of the crystal structure and atoms involved. However, to maximize the ionic conductivity, several methods 
        can be done. Firstly, operating at higher temperatures can significantly decrease these ohmic losses. 
        Substitutional doping methods to further refine the crystal structure and control defect concentrations 
        can also play a significant role in increasing the conductivity. Another way to decrease ohmic resistance 
        is to decrease the thickness of the electrolyte layer.

        An ionic specific resistance of the electrolyte as a function of temperature can be described by the following relationship:
        ''')

        st.markdown(r'''
        The general expression for theohmic loss is given by Ohm’s law in a finite form:
        ''')

        st.latex(r''' E_{\text{Ohm}}=R_{\text{Ohm}} \cdot I''')

        st.markdown(r'''
        where the resistance $R_{\text{Ohm}}$ depended from internal parameters of fuel cell materials as a current conductor, specsificaly $\sigma$ (average area conductivity) and
        $d$ (elrolyte thickness).
        ''')


        st.markdown(r'''
        In order to find out how $E_{\text{Ohm}}$  and $R_{\text{Ohm}}$ depends on $\sigma$ and  $\delta$ complete the following tasks:
        ''')


        st.latex(r''' r_1=\frac{\delta}{\sigma}''')

        st.markdown('''
        where:  
        >$\delta$ - electrolyte thickness  
        >$\sigma$ - ionic conductivity  

        The ionic conductivity of the solid oxide is defined as follows:
        ''')


        st.latex(r''' \sigma = \sigma_0 \cdot e^{\frac{-E}{R \cdot T}}''')

        st.markdown('''
        where:  
        >$\sigma_{0}$ and $E$ – factors depended on electrolyte materials  
        >$T$ – electrolyte temperature  
        >$R$ – ideal gas constant  
        ''')


        st.markdown('''
        Complete the following tasks: 
        ''')

    with st.container():

        st.markdown(task_title("Задание 2.1","Electrolyte conductivity and its role in Ohmic losses"), unsafe_allow_html=True)

        st.markdown(task_subtitle("Работа с моделью"), unsafe_allow_html=True)

        st.markdown(r'''
        * Using the sidebars at the `Model parameters` panel set parameter $\sigma$  to 2.5 S/m. 
        * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
        * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.
        ''')



        col1, col2 = st.columns([2,1])
        chk_vi5 = col2.checkbox(td[lang_selected][14], key="chk_vi5")
        if chk_vi5 and btn_runsimulation:
            st.session_state["vi5"] = [
                temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
        st.bokeh_chart(plt_vi5.makePlotVI(st.session_state["vi5"]))
        col1.markdown(plt_statetitle("vi5"), unsafe_allow_html=True)


        # with st.expander("График"):
        #     col1, col2 = st.columns([3,1])
        #     chk_cond1 = col2.checkbox("Allow to refresh", key="chk_cond1")
        #     if chk_cond1 and btn_runsimulation:
        #         st.session_state["svi1"] = [
        #             temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
        #     st.bokeh_chart(makePlotVI(st.session_state["svi1"]))
        #     col1.write("Saved state: T = {0:.0f} K, σ = {2:.2f}⋅ 10⁻⁶ S/m, d = {3:.2f} µm, jₘ = {4:.2f} mA/cm²".format(*st.session_state["svi1"]))

        st.markdown('''
        * Using the sidebars at the `Model parameters` panel set parameter $\sigma$  to 0.2 S/m. 
        * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
        * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.
        ''')


        col1, col2 = st.columns([2,1])
        chk_vi6 = col2.checkbox(td[lang_selected][14], key="chk_vi6")
        if chk_vi6 and btn_runsimulation:
            st.session_state["vi6"] = [
                temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
        st.bokeh_chart(plt_vi6.makePlotVI(st.session_state["vi6"]))
        col1.markdown(plt_statetitle("vi6"), unsafe_allow_html=True)

        # with st.expander("График"):
        #     col1, col2 = st.columns([3,1])
        #     chk_cond1 = col2.checkbox("Allow to refresh", key="chk_cond2")
        #     if chk_cond1 and btn_runsimulation:
        #         st.session_state["svi2"] = [
        #             temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
        #     st.bokeh_chart(makePlotVI(st.session_state["svi2"]))
        #     col1.write("Saved state: T = {0:.0f} K, σ = {2:.2f}⋅ 10⁻⁶ S/m, d = {3:.2f} µm, jₘ = {4:.2f} mA/cm²".format(*st.session_state["svi2"]))



        st.markdown(task_subtitle("Анализ результатов"), unsafe_allow_html=True)

        st.markdown('''
        Compare the obtained Voltage-Current curves  and answer the questions:
         
        * How the increasing of conductivity affects on cell performance?

        ''')

        
    with st.container():    

        st.markdown(task_title("Задание 2.2","Electrolyte thickness and its role in  Ohmic losses"), unsafe_allow_html=True)

        st.markdown(task_subtitle("Работа с моделью"), unsafe_allow_html=True)

        st.markdown(r'''
        * Using the sidebars at the `Model parameters` panel set parameter $d$  to 30 µm. 
        * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
        * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.
        ''')

        st.markdown(
            '''
            <style>
            .reportview-container {
                background: url("url_goes_here")
            }
            .sidebar .sidebar-content {
                background: url("url_goes_here")
            }
            </style>
            ''',
            unsafe_allow_html=True
        )


        col1, col2 = st.columns([2,1])
        chk_vi7 = col2.checkbox(td[lang_selected][14], key="chk_vi7")
        if chk_vi7 and btn_runsimulation:
            st.session_state["vi7"] = [
                temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
        st.bokeh_chart(plt_vi7.makePlotVI(st.session_state["vi7"]))
        col1.markdown(plt_statetitle("vi7"), unsafe_allow_html=True)

        # with st.expander("График"):
        #     col1, col2 = st.columns([3,1])
        #     chk_cond1 = col2.checkbox("Allow to refresh", key="chk_cond3")
        #     if chk_cond1 and btn_runsimulation:
        #         st.session_state["svi3"] = [
        #             temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
        #     st.bokeh_chart(makePlotVI(st.session_state["svi3"]))
        #     col1.write("Saved state: T = {0:.0f} K, σ = {2:.2f}⋅ 10⁻⁶ S/m, d = {3:.2f} µm, jₘ = {4:.2f} mA/cm²".format(*st.session_state["svi3"]))



        st.markdown(r'''
        * Using the sidebars at the `Model parameters` panel set parameter  $d$ to 150 µm.
        * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
        * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.
        ''')

        col1, col2 = st.columns([2,1])
        chk_vi8 = col2.checkbox(td[lang_selected][14], key="chk_vi8")
        if chk_vi8 and btn_runsimulation:
            st.session_state["vi8"] = [
                temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
        st.bokeh_chart(plt_vi8.makePlotVI(st.session_state["vi8"]))
        col1.markdown(plt_statetitle("vi8"), unsafe_allow_html=True)


        # with st.expander("График"):
        #     col1, col2 = st.columns([3,1])
        #     chk_cond1 = col2.checkbox("Allow to refresh", key="chk_cond4")
        #     if chk_cond1 and btn_runsimulation:
        #         st.session_state["svi4"] = [
        #             temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
        #     st.bokeh_chart(makePlotVI(st.session_state["svi4"]))
        #     col1.write("Saved state: T = {0:.0f} K, σ = {2:.2f}⋅ 10⁻⁶ S/m, d = {3:.2f} µm, jₘ = {4:.2f} mA/cm²".format(*st.session_state["svi4"]))

        st.markdown(task_subtitle("Анализ результатов"), unsafe_allow_html=True)

        st.markdown('''
        Compare the obtained Voltage-Current curves  and answer the questions:
         
        * How the increasing of electrolyte thickness affects on cell performance?
        * What is the relationship between $d$  and $\sigma$ at fixed $V$?


        ''')


if tab_selected == works[2] and lang_selected == langs[1]:

    cflag = True 

    if cflag:
        st.markdown('''
        Under construction...
        ''')
    else:

        with st.container():

            st.header("Ohmic losses")

            st.markdown(r'''
            Ohmic losses ($E_{\text{Ohm}}$) in an SOFC result from ionic conductivity through the electrolyte and electrical resistance 
            offered to the flow of electrons in the external electrical circuit. This is inherently a materials property 
            of the crystal structure and atoms involved. However, to maximize the ionic conductivity, several methods 
            can be done. Firstly, operating at higher temperatures can significantly decrease these ohmic losses. 
            Substitutional doping methods to further refine the crystal structure and control defect concentrations 
            can also play a significant role in increasing the conductivity. Another way to decrease ohmic resistance 
            is to decrease the thickness of the electrolyte layer.

            An ionic specific resistance of the electrolyte as a function of temperature can be described by the following relationship:
            ''')

            st.markdown(r'''
            The general expression for theohmic loss is given by Ohm’s law in a finite form:
            ''')

            st.latex(r''' E_{\text{Ohm}}=R_{\text{Ohm}} \cdot I''')

            st.markdown(r'''
            where the resistance $R_{\text{Ohm}}$ depended from internal parameters of fuel cell materials as a current conductor, specsificaly $\sigma$ (average area conductivity) and
            $d$ (elrolyte thickness).
            ''')


            st.markdown(r'''
            In order to find out how $E_{\text{Ohm}}$  and $R_{\text{Ohm}}$ depends on $\sigma$ and  $\delta$ complete the following tasks:
            ''')


            st.latex(r''' r_1=\frac{\delta}{\sigma}''')

            st.markdown('''
            where:  
            >$\delta$ - electrolyte thickness  
            >$\sigma$ - ionic conductivity  

            The ionic conductivity of the solid oxide is defined as follows:
            ''')


            st.latex(r''' \sigma = \sigma_0 \cdot e^{\frac{-E}{R \cdot T}}''')

            st.markdown('''
            where:  
            >$\sigma_{0}$ and $E$ – factors depended on electrolyte materials  
            >$T$ – electrolyte temperature  
            >$R$ – ideal gas constant  
            ''')


            st.markdown('''
            Complete the following tasks: 
            ''')


            with st.expander("Task 1"):

                st.subheader("1. Electrolyte conductivity and its role in Ohmic losses")
                st.markdown(r'''
                * Using the sidebars at the `Model parameters` panel set parameter $\sigma$ to 2.5 S/m. 
                * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
                * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.

                ---

                ''')

                col1, col2 = st.columns([3,1])
                chk_cond1 = col2.checkbox("Allow to refresh", key="chk_cond1")
                if chk_cond1 and btn_runsimulation:
                    st.session_state["svi1"] = [
                        temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
                st.bokeh_chart(makePlotVI(st.session_state["svi1"]))
                col1.write("Saved state: T = {0:.0f} K, σ = {2:.2f}⋅ 10⁻⁶ S/m, d = {3:.2f} µm, jₘ = {4:.2f} mA/cm²".format(*st.session_state["svi1"]))

                st.markdown('''
                * Using the sidebars at the `Model parameters` panel set parameter $\sigma$  to 0.2 S/m. 
                * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
                * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.

                ---

                ''')

                
                col1, col2 = st.columns([3,1])
                chk_cond1 = col2.checkbox("Allow to refresh", key="chk_cond2")
                if chk_cond1 and btn_runsimulation:
                    st.session_state["svi2"] = [
                        temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
                st.bokeh_chart(makePlotVI(st.session_state["svi2"]))
                col1.write("Saved state: T = {0:.0f} K, σ = {2:.2f}⋅ 10⁻⁶ S/m, d = {3:.2f} µm, jₘ = {4:.2f} mA/cm²".format(*st.session_state["svi2"]))

                st.markdown('''    

                ---

                ''')



                st.subheader("2. Results analysis")

                st.markdown('''
                Compare the obtained Voltage-Current curves  and answer the questions:
                 
                * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
                * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.

                ''')

            with st.expander("Task 2"):

                st.subheader("1. Electrolyte thickness and its role in  Ohmic losses")
                st.markdown(r'''
                * Using the sidebars at the `Model parameters` panel set parameter $d$ to 30 µm. 
                * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
                * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.

                ---

                ''')

                st.markdown(
                    '''
                    <style>
                    .reportview-container {
                        background: url("url_goes_here")
                    }
                    .sidebar .sidebar-content {
                        background: url("url_goes_here")
                    }
                    </style>
                    ''',
                    unsafe_allow_html=True
                )

                with st.container():

                    col1, col2 = st.columns([3,1])
                    chk_cond1 = col2.checkbox("Allow to refresh", key="chk_cond3")
                    if chk_cond1 and btn_runsimulation:
                        st.session_state["svi3"] = [
                            temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
                    st.bokeh_chart(makePlotVI(st.session_state["svi3"]))
                    col1.write("Saved state: T = {0:.0f} K, σ = {2:.2f}⋅ 10⁻⁶ S/m, d = {3:.2f} µm, jₘ = {4:.2f} mA/cm²".format(*st.session_state["svi3"]))

                st.markdown('''    

                ---

                ''')


                st.markdown(r'''
                * Using the sidebars at the `Model parameters` panel set parameter  $d$ to 30 µm.
                * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
                * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.

                ---

                ''')

                
                col1, col2 = st.columns([3,1])
                chk_cond1 = col2.checkbox("Allow to refresh", key="chk_cond4")
                if chk_cond1 and btn_runsimulation:
                    st.session_state["svi4"] = [
                        temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
                st.bokeh_chart(makePlotVI(st.session_state["svi4"]))
                col1.write("Saved state: T = {0:.0f} K, σ = {2:.2f}⋅ 10⁻⁶ S/m, d = {3:.2f} µm, jₘ = {4:.2f} mA/cm²".format(*st.session_state["svi4"]))

                st.markdown('''    

                ---

                ''')



                st.subheader("2. Results analysis")

                st.markdown('''
                Compare the obtained Voltage-Current curves  and answer the questions:
                 
                * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
                * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.

                ''')


if tab_selected == works[3] and lang_selected == langs[0]:

    with st.container():

        st.header("Mass transport losses")
 
        st.markdown('''
        <div style="text-align: justify;">        
              To produce electricity, a fuel cell must be continually supplied with fuel and oxidant. 
        At the same time, products must be continuously removed so as to avoid “strangling” the cell. 
        The process of supplying reactants and removing products is termed fuel cell mass transport. 
        As you will learn, this seemingly simple task can turn out to be quite complicated.

        Why are we so interested in fuel cell mass transport? The answer is because poor mass transport leads to significant fuel cell performance losses. 
        To understand why poor mass transport can lead to a performance loss, remember that fuel cell performance 
        is determined by the reactant and product concentrations within the catalyst layer, not at the fuel cell inlet. 
        Thus, reactant depletion (or product accumulation) within the catalyst layer will adversely affect performance. 
        This loss in performance is called a fuel cell *concentration* loss.

        Concentration polarization losses are sometimes expressed as a function of the limiting current, $j_m$,
        which is usually taken as a measure of the maximum rate atwhich a reactant can be supplied to an electrode:
        </div>
        ''', unsafe_allow_html=True)

        st.latex(
            r'''
            E_{con} = \frac{RT}{2F}\cdot\ln\left(\frac{j_m}{j_m-j}\right)
            '''
        )


        st.markdown('''
        <div style="text-align: justify;">        
        To understand the influnce of the limiting current on the performance of the fuel cell complete the tollowing task.
        <p></p></div>
        ''', unsafe_allow_html=True)


    with st.container():

        st.markdown(task_title("Задание 3.1","Test a low and high limiting current"), unsafe_allow_html=True)

        st.markdown(task_subtitle("Работа с моделью"), unsafe_allow_html=True)

    


        st.markdown('''
        * Using the sidebars at the `Model parameters` panel change cell limiting current to 60 mA/cm². 
        * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
        * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.
        ''')


        col1, col2 = st.columns([2,1])
        chk_pi1 = col2.checkbox(td[lang_selected][14], key="chk_pi1")
        if chk_pi1 and btn_runsimulation:
            st.session_state["pi1"] = [
                temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
        st.bokeh_chart(plt_pi1.makePlotPI(st.session_state["pi1"]))
        col1.markdown(plt_statetitle("pi1"), unsafe_allow_html=True)



        # col1, col2 = st.columns([3,1])
        # chk_temp2 = col2.checkbox("Allow to refresh", key="chk_mass1")
        # if chk_temp2 and btn_runsimulation:
        #     st.session_state["mvi1"] = [
        #         temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
        # st.bokeh_chart(makePlotPI(st.session_state["mvi1"]))
        # col1.write("Saved state: T = {0:.0f} K, σ = {2:.2f}⋅ 10⁻⁶ S/m, d = {3:.2f} µm, jₘ = {4:.2f} mA/cm²".format(*st.session_state["mvi1"]))


        st.markdown('''
        * Using the sidebars at the `Model parameters` panel change cell limiting current to 100 mA/cm².  
        * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
        * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.
        ''')

        col1, col2 = st.columns([2,1])
        chk_pi2 = col2.checkbox(td[lang_selected][14], key="chk_pi2")
        if chk_pi2 and btn_runsimulation:
            st.session_state["pi2"] = [
                temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
        st.bokeh_chart(plt_pi2.makePlotPI(st.session_state["pi2"]))
        col1.markdown(plt_statetitle("pi2"), unsafe_allow_html=True)
        

        # col1, col2 = st.columns([3,1])
        # chk_temp2 = col2.checkbox("Allow to refresh", key="chk_mass2")
        # if chk_temp2 and btn_runsimulation:
        #     st.session_state["mvi2"] = [
        #         temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
        # st.bokeh_chart(makePlotPI(st.session_state["mvi2"]))
        # col1.write("Saved state: T = {0:.0f} K, σ = {2:.2f}⋅ 10⁻⁶ S/m, d = {3:.2f} µm, jₘ = {4:.2f} mA/cm²".format(*st.session_state["mvi2"]))



        st.markdown(task_subtitle("Анализ результатов"), unsafe_allow_html=True)

        st.markdown('''
        Compare the obtained Power-Current curves  and answer the questions:
         
        * How the limiting current density affects on cell performance?


        ''')


if tab_selected == works[3] and lang_selected == langs[1]:

    cflag = True 

    if cflag:
        st.markdown('''
        Under construction...
        ''')
    else:


        with st.container():

            st.header("Mass transport losses")
 
            st.markdown('''
            <div style="text-align: justify;">        
                  To produce electricity, a fuel cell must be continually supplied with fuel and oxidant. 
            At the same time, products must be continuously removed so as to avoid “strangling” the cell. 
            The process of supplying reactants and removing products is termed fuel cell mass transport. 
            As you will learn, this seemingly simple task can turn out to be quite complicated.

            Why are we so interested in fuel cell mass transport? The answer is because poor mass transport leads to significant fuel cell performance losses. 
            To understand why poor mass transport can lead to a performance loss, remember that fuel cell performance 
            is determined by the reactant and product concentrations within the catalyst layer, not at the fuel cell inlet. 
            Thus, reactant depletion (or product accumulation) within the catalyst layer will adversely affect performance. 
            This loss in performance is called a fuel cell *concentration* loss.

            Concentration polarization losses are sometimes expressed as a function of the limiting current, $j_m$,
            which is usually taken as a measure of the maximum rate atwhich a reactant can be supplied to an electrode:
            </div>
            ''', unsafe_allow_html=True)

            st.latex(
                r'''
                E_{con} = \frac{RT}{2F}\cdot\ln\left(\frac{j_m}{j_m-j}\right)
                '''
            )


            st.markdown('''
            <div style="text-align: justify;">        
            To understand the influnce of the limiting current on the performance of the fuel cell complete the tollowing task.
            <p></p></div>
            ''', unsafe_allow_html=True)


            with st.expander("Task 1"):

                st.subheader("1. Test a low and high limiting current")

                st.markdown('''
                * Using the sidebars at the `Model parameters` panel change cell temperature to 750 K. 
                * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
                * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.

                ---

                ''')

                col1, col2 = st.columns([3,1])
                chk_temp2 = col2.checkbox("Allow to refresh", key="chk_mass1")
                if chk_temp2 and btn_runsimulation:
                    st.session_state["mvi1"] = [
                        temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
                st.bokeh_chart(makePlotPI(st.session_state["mvi1"]))
                col1.write("Saved state: T = {0:.0f} K, σ = {2:.2f}⋅ 10⁻⁶ S/m, d = {3:.2f} µm, jₘ = {4:.2f} mA/cm²".format(*st.session_state["mvi1"]))

                st.markdown('''    

                ---

                ''')

                st.markdown('''
                * Using the sidebars at the `Model parameters` panel change cell temperature to 750 K. 
                * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
                * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.

                ---

                ''')

                col1, col2 = st.columns([3,1])
                chk_temp2 = col2.checkbox("Allow to refresh", key="chk_mass2")
                if chk_temp2 and btn_runsimulation:
                    st.session_state["mvi2"] = [
                        temperature, pressure, sigma, ethick, jm, H2ac, H2Oac, O2cc]
                st.bokeh_chart(makePlotPI(st.session_state["mvi2"]))
                col1.write("Saved state: T = {0:.0f} K, σ = {2:.2f}⋅ 10⁻⁶ S/m, d = {3:.2f} µm, jₘ = {4:.2f} mA/cm²".format(*st.session_state["mvi2"]))

                
                st.markdown('''    

                ---

                ''')



                st.subheader("2. Results analysis")

                st.markdown('''
                Compare the obtained Voltage-Current curves  and answer the questions:
                 
                * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
                * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.

                ''')


if "ni_fuel_gas" not in st.session_state.keys():
    st.session_state['ni_fuel_gas'] = np.round(10.0,1)

if "ni_eff_value" not in st.session_state.keys():
    st.session_state['ni_eff_value'] = np.round(65.0,1)
if "ni_gen_power" not in st.session_state.keys():
    st.session_state['ni_gen_power'] = np.round(216.0,1)

if "ni_gen_adduct" not in st.session_state.keys():
    st.session_state['ni_gen_adduct'] = np.round(89.0,1)



def _update_on_change_fuel_gas():
    # st.session_state['ni_fuel_gas'] = st.session_state['ni_fuel_gas']
    # st.session_state['ni_eff_value'] = eff_value
    st.session_state['ni_gen_power'] = st.session_state['ni_fuel_gas']*st.session_state['ni_eff_value']/100.0*33.0
    st.session_state['ni_gen_adduct'] = st.session_state['ni_fuel_gas']*8.9360119
def _update_on_change_eff_value():
    # st.session_state['ni_fuel_gas'] = st.session_state['ni_fuel_gas']
    # st.session_state['ni_eff_value'] = eff_value
    # st.session_state['ni_gen_power'] = st.session_state['ni_fuel_gas']*st.session_state['ni_eff_value']/100.0*33.0
    # st.session_state['ni_gen_adduct'] = st.session_state['ni_fuel_gas']*8.9360119
    _update_on_change_fuel_gas()
def _update_on_change_gen_power():
    st.session_state['ni_fuel_gas'] = st.session_state['ni_gen_power']/st.session_state['ni_eff_value']*100.0/33.326
    # st.session_state['ni_eff_value'] = eff_value
    # st.session_state['ni_gen_power'] = st.session_state['ni_fuel_gas']*st.session_state['ni_eff_value']/100.0*33.0
    st.session_state['ni_gen_adduct'] = st.session_state['ni_fuel_gas']*8.9360119
def _update_on_change_gen_adduct():
    st.session_state['ni_fuel_gas'] = st.session_state['ni_gen_adduct']/8.9360119
    # st.session_state['ni_eff_value'] = eff_value
    st.session_state['ni_gen_power'] = st.session_state['ni_fuel_gas']*st.session_state['ni_eff_value']/100.0*33.326
    # st.session_state['ni_gen_adduct'] = st.session_state['ni_fuel_gas']*8.9360119




if tab_selected == works[4] and lang_selected == langs[0]:
    
    with st.container():
        st.header("Fuel Cell Calculator")
        st.markdown('''
        <div style="text-align: justify;">        
        Calculate the energy and the amount of water a fuel cell produces when burning hydrogen. 
        A fuel cell is an energy converter, which transforms chemical binding energy into electricity. 
        Often hydrogen, H₂, is used as a fuel and transformed into water, H₂O. This calculator refers to the use with hydrogen. 
        Fuel cells can also use other materials, like methane.
        <p></p></div>
        ''', unsafe_allow_html=True)



        # fuel_gas, eff_value, gen_power, gen_adduct = 0.0, 65.0, 0.0, 0.0




        col1, col2 = st.columns([1,1]) 
        fuel_gas = col1.number_input(td[lang_selected][20], value=np.round(st.session_state['ni_fuel_gas'],1), format="%.1f", min_value=0.0, step=0.1, key="ni_fuel_gas", on_change=_update_on_change_fuel_gas)
        eff_value = col2.number_input(td[lang_selected][21], value=np.round(st.session_state['ni_eff_value'],1), format="%.0f",  min_value=0.1, max_value=100.0, step=1.0, key="ni_eff_value", on_change=_update_on_change_eff_value)
        col1, col2 = st.columns([1,1])
        gen_power = col1.number_input(td[lang_selected][22], value=np.round(st.session_state['ni_gen_power'],1), format="%.1f", min_value=0.0, step=0.1, key="ni_gen_power", on_change=_update_on_change_gen_power)
        gen_adduct = col2.number_input(td[lang_selected][23], value=np.round(st.session_state['ni_gen_adduct'],1), format="%.1f", min_value=0.0, step=0.1, key="ni_gen_adduct", on_change=_update_on_change_gen_adduct)


        # st.session_state['fuel_gas'] = fuel_gas
        # st.session_state['eff_value'] = eff_value
        # st.session_state['gen_power'] = gen_power
        # st.session_state['gen_adduct'] = gen_adduct




        st.markdown(task_title("Задание 4.1","Расчет характеристик энергоустановки"), unsafe_allow_html=True)


        st.markdown('''
            * Используя представленный выше калькулятор, посчитайте, какое количество топлива (газа водорода) понадобится для снабжения электроэнергией двухкомнатной квартры в которой проживает семья из трех человека? Средняя потребность одного человека в электроэнергии составляет около 63 кВтч в месяц. 
            * Какова будет стоимость месяца работы такой энергоустановки (средняя цена 1 кг водорода 800-1000 руб.)?
            * Можно ли оценить, какое количество тепла будет генерировать такая энергоустановка?
            ''')
        


if tab_selected == works[4] and lang_selected == langs[1]:

    cflag = True 

    if cflag:
        st.markdown('''
        Under construction...
        ''')
    else:
    
        with st.container():
            st.header("Fuel Cell Calculator")
            st.markdown('''
            <div style="text-align: justify;">        
            Calculate the energy and the amount of water a fuel cell produces when burning hydrogen. 
            A fuel cell is an energy converter, which transforms chemical binding energy into electricity. 
            Often hydrogen, H₂, is used as a fuel and transformed into water, H₂O. This calculator refers to the use with hydrogen. 
            Fuel cells can also use other materials, like methane.
            <p></p></div>
            ''', unsafe_allow_html=True)


            with st.expander("Task 1"):

                st.markdown('''
                    * Using the sidebars at the `Model parameters` panel change cell temperature to 750 K. 
                    * Check `Allow to refresh` flag and press `Run simulation` button under the panel.
                    * Plot in figure below should update. Uncheck `Allow to refresh` flag to forbid any chages of the plot.
                    ''')
            

                col1, col2 = st.columns([1,1]) 
                fuel_gas = col1.number_input("Total consumed H₂, [kg]", value=np.round(st.session_state['ni_fuel_gas'],1), format="%.1f", min_value=0.0, step=0.1, key="ni_fuel_gas", on_change=_update_on_change_fuel_gas)
                eff_value = col2.number_input("Efficiency factor, [%]", value=np.round(st.session_state['ni_eff_value'],1), format="%.1f",  min_value=0.1, max_value=100.0, step=0.1, key="ni_eff_value", on_change=_update_on_change_eff_value)
                col1, col2 = st.columns([1,1])
                gen_power = col1.number_input("Generated power, [kWh]", value=np.round(st.session_state['ni_gen_power'],1), format="%.1f", min_value=0.0, step=0.1, key="ni_gen_power", on_change=_update_on_change_gen_power)
                gen_adduct = col2.number_input("Generated H₂O, [kg]", value=np.round(st.session_state['ni_gen_adduct'],1), format="%.1f", min_value=0.0, step=0.1, key="ni_gen_adduct", on_change=_update_on_change_gen_adduct)


                # st.session_state['fuel_gas'] = fuel_gas
                # st.session_state['eff_value'] = eff_value
                # st.session_state['gen_power'] = gen_power
                # st.session_state['gen_adduct'] = gen_adduct
    

 



# if tab_selected == works[6]:
#     introduction_container = st.container()


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
