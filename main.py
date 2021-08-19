from pathlib import WindowsPath
import numpy as np
from requests.api import head
import streamlit as st
import requests
import pandas as pd
from bokeh.plotting import figure

#st.set_page_config(page_title=None, page_icon=None, layout='wide', initial_sidebar_state='auto')

API_URL = 'http://178.154.215.108/tote/ivc'
works = ['Introduction', 
        'Work #1: Understanding IVC', 
        'Work #2: The second work', 
        'Work #3: The third work',
        'Futher investigation']

# ==============================
# Sidebar
# ==============================

st.sidebar.image('./Images/logo.png')
st.sidebar.markdown('***')
st.sidebar.header('SOFC cloud model')
tab_selected = st.sidebar.selectbox('Select the work:', works)

if tab_selected != works[0]:
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
    
# ==============================

def get_ivc():
    payload = {
            'temperature': temperature,
            'sigma': sigma,
            'ethick': ethick * pow(10, -6)
        }
    result = requests.post(API_URL, data=payload)
    return result.json()

# ==============================
# Main page
# ==============================

header_container = st.beta_container()
for item in works:
    if tab_selected == item:
        with header_container:
            st.title(item)
            '''***'''

if tab_selected == works[0]:
    
    introduction_container = st.beta_container()
    scheme_container = st.beta_container()
    basics_container = st.beta_container()
    theory_container = st.beta_container()
    whatnext_container = st.beta_container()

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
        st.latex(r''' cathode: \frac{1}{2}O_2+2e^- \rightarrow O^{2-}''')
        with st.beta_expander('Anode description'):
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
        st.latex(r''' anode: H_2+O^{2-} \rightarrow H_2O+2e^-''')
        with st.beta_expander('Cathode description'):
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

        with st.beta_expander('Electrolyte description'):
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
        with st.beta_expander('Interconnect'):
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
        st.latex(r'''E_{SOFC}=\frac{E_{max}-i_{max} \cdot \eta_f\cdot r_1}{\frac{r_1}{r_2}\cdot (1-\eta_f)+1}''')
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
        with st.beta_expander('Ohmic polarization'):
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
        with st.beta_expander('Concentration polarization'):
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
        with st.beta_expander('Activation polarization'):
            '''
            The activation polarization is the result of the kinetics involved with the electrochemical reactions. 
            Each reaction has a certain activation barrier that must be overcome in order to proceed and this barrier 
            leads to the polarization. The activation barrier is the result of many complex electrochemical reaction steps 
            where typically the rate limiting step is responsible for the polarization. The polarization equation shown below 
            is found by solving the [Butler–Volmer](https://en.wikipedia.org/wiki/Butler–Volmer_equation) equation in the 
            high current density regime (where the cell typically operates), and can be used to estimate the activation polarization:
            '''
            st.latex(r''' \eta_{act}=\frac{RT}{\beta zF} \cdot ln \left( \frac{i}{i_{0}} \right)''')
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

if tab_selected == works[1]:
    
    introduction_container = st.beta_container()


    with introduction_container:
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
        of $1.229 V$ at standard conditions $(T = 273.15 K, p = 1 atm)$. 
        With the standard electrode potential $E_0$, universal gas constant $R$, 
        temperature $T$, Faraday’s constant $F$, molar concentration $x$ and pressure $p$, 
        the OCV is given by the following:
        '''
        st.latex(r'''E(p,T)=E^0(p^0,T)-\frac{RT}{2F}\ln \left( \frac{x_{H_2O}}{x_{H_2} \cdot \sqrt{x_{O_2}}} \right) + \frac{RT}{4F}\ln{p}''')
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

    if st.sidebar.button('Calculate'):
        result = get_ivc()
        p = figure()
        p.line(result['i'],result['Eload'])
        st.bokeh_chart(p, use_container_width=True)


if tab_selected == works[2]:
    introduction_container = st.beta_container()

if tab_selected == works[3]:
    introduction_container = st.beta_container()

if tab_selected == works[4]:
    introduction_container = st.beta_container()


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
