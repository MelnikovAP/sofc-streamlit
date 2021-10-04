#+---------------------+---------------+-----------------------------+
#|   Параметры модели  |  Обозначение  | Область допустимых значений |
#|                     |               |        и форма записи       |
#+=====================+===============+=============================+
#| Температура         |       T       |        1000+/-20% [K]       |
#+---------------------+---------------+-----------------------------+
#| Давление газовой    |       P       |        101325.0 [Pa]        |
#| смеси               |               |                             |
#+---------------------+---------------+-----------------------------+
#| Состав анодной      |  anode_gas_X  |     'H2:0.97, H2O:0.03'     |
#| газовой смеси       |               |                             |
#+---------------------+---------------+-----------------------------+
#| Состав катодной     | cathode_gas_X |     'O2:1.0, H2O:0.001'     |
#| газовой смеси       |               |                             |
#+---------------------+---------------+-----------------------------+
#| Проводимость        |     sigma     |       2.0+/-10% [См/м]      |
#| электролита         |               |                             |
#+---------------------+---------------+-----------------------------+
#| Толщина электролита |     ethick    |       5.0e-5+/-10% [м]      |
#+---------------------+---------------+-----------------------------+
import cantera as ct
import numpy as np
import logging
from scipy.optimize import fsolve

tsofclogger = logging.getLogger("tsofc")
ct.add_module_directory()


class TSOFC(object):
    def __init__(self, **kwargs):
        # Настраиваемые параметры
        self.T = kwargs.get("T", 1073.15)  # T в кельвинах 
        self.P = kwargs.get("P", ct.one_atm)  # 1 атм. в Па =  101325.0 Па
        self.anode_gas_X = kwargs.get("anode_gas_X", 'H2:0.97, H2O:0.03') # составы газовых смесей в долях
        self.cathode_gas_X = kwargs.get("cathode_gas_X",'O2:1.0, H2O:0.001')
        self.sigma = kwargs.get("sigma", 2.0)  # проводимость электролита [См / м]
        self.ethick = kwargs.get("ethick", 5.0e-5)  # толщина электролита [м]
        
        # эмпирические параметры для учета концетрационных потерь
        # self.m  = kwargs.get("m", 3.0e-5)
        self.jm = kwargs.get("jm", 100)
        self.alpha = 0.5

        # Внутренние параметры
        self.tss = 50.0
        self.TPB_length_per_area = 1.0e7  # размер трехфазной области на единицу площади [1/m]

        self.gas_a, self.anode_bulk, self.oxide_a = ct.import_phases(
            'sofc.yaml', ['gas', 'metal', 'oxide_bulk'])
        self.anode_surf = ct.Interface('sofc.yaml', 'metal_surface', [self.gas_a])
        self.oxide_surf_a = ct.Interface('sofc.yaml', 'oxide_surface', [self.gas_a, self.oxide_a])
        self.tpb_a = ct.Interface('sofc.yaml', 'tpb', [self.anode_bulk, self.anode_surf, self.oxide_surf_a])
        self.anode_surf.name = 'anode surface'
        self.oxide_surf_a.name = 'anode-side oxide surface'

        self.gas_c, self.cathode_bulk, self.oxide_c = ct.import_phases(
            'sofc.yaml', ['gas', 'metal', 'oxide_bulk'])
        self.cathode_surf = ct.Interface('sofc.yaml', 'metal_surface', [self.gas_c])
        self.oxide_surf_c = ct.Interface('sofc.yaml', 'oxide_surface', [self.gas_c, self.oxide_c])
        self.tpb_c = ct.Interface('sofc.yaml', 'tpb', [self.cathode_bulk, self.cathode_surf, self.oxide_surf_c])
        self.cathode_surf.name = 'cathode surface'
        self.oxide_surf_c.name = 'cathode-side oxide surface'

        self.gas_a.TPX = self.T, self.P, self.anode_gas_X
        self.gas_a.equilibrate('TP')
        self.gas_c.TPX = self.T, self.P, self.cathode_gas_X
        self.gas_c.equilibrate('TP')

        phases = [
            self.anode_bulk, self.anode_surf, self.oxide_surf_a, self.oxide_a, 
            self.cathode_bulk, self.cathode_surf, self.oxide_surf_c, self.oxide_c, 
            self.tpb_a, self.tpb_c
        ]
        for p in phases:
            p.TP = self.T, self.P

        for s in [self.anode_surf, self.oxide_surf_a, self.cathode_surf, self.oxide_surf_c]:
            s.advance_coverages(self.tss)
            show_coverages(s)

        self.Ea0 = Solver(self.anode_curr, xstart=-0.51)
        self.Ec0 = Solver(self.cathode_curr, xstart=0.51)

        tsofclogger.info('\nocv from zero current is: {0:.3e}'.format(float((self.Ec0 - self.Ea0).flat[0])))
        tsofclogger.info('OCV from thermo equil is: {0:.3e}'.format(float(equil_OCV(self.gas_a, self.gas_c).flat[0])))
        tsofclogger.info('Ea0 = = {0:.3e}'.format(float(self.Ea0.flat[0])))
        tsofclogger.info('Ec0 = {0:.3e}'.format(float(self.Ec0.flat[0])))


    def anode_curr(self, E):
        self.anode_bulk.electric_potential = E
        w = self.tpb_a.net_production_rates
        electron_index = self.tpb_a.kinetics_species_index('electron')
        return ct.faraday * w[electron_index] * self.TPB_length_per_area

    def cathode_curr(self, E):
        self.cathode_bulk.electric_potential = E + self.oxide_c.electric_potential
        w = self.tpb_c.net_production_rates
        electron_index = self.tpb_c.kinetics_species_index('electron')
        return -ct.faraday * w[electron_index] * self.TPB_length_per_area

    def getPolarizationCurve(self):   
        Ea_min = self.Ea0 - 0.25
        Ea_max = self.Ea0 + 0.25
        output_data = []
        #for n in range(350):
        for n in range(400):
            Ea = Ea_min + 0.005*n
            self.anode_bulk.electric_potential = Ea
            curr = self.anode_curr(Ea)
            delta_V = curr * self.ethick / self.sigma
            phi_oxide_c = -delta_V
            self.oxide_c.electric_potential = phi_oxide_c
            self.oxide_surf_c.electric_potential = phi_oxide_c
            Ec = Solver(self.cathode_curr, xstart=self.Ec0+0.1, C=curr)
            # E_mass = self.m * np.exp(self.n * curr)

            eta_con = (1.0 + 1.0 / self.alpha) * ct.gas_constant * self.T / \
                (2.0 * ct.faraday) * np.log(self.jm/(self.jm-0.1*curr))

            self.cathode_bulk.electric_potential = phi_oxide_c + Ec
            output_data.append([float(x) for x in [0.1*curr, Ea - self.Ea0, Ec - self.Ec0, delta_V,
                                self.cathode_bulk.electric_potential -
                                self.anode_bulk.electric_potential, eta_con]])
        output_data = np.array(output_data)                            
        return output_data[(output_data[:,0]>0)&(output_data[:,-1]>0)] # Выходные данные в формате 'i (mA/cm2)', 'eta_a', 'eta_c', 'eta_ohmic', 'Eload'


def show_coverages(s):
     tsofclogger.info('\n{0:s}\n'.format(s.name))
     cov = s.coverages
     names = s.species_names
     for n in range(s.n_species):
         tsofclogger.info('{0:16s}  {1:13.4g}'.format(names[n], cov[n]))


def equil_OCV(gas1, gas2):
    return (-ct.gas_constant * gas1.T *
            np.log(gas1['O2'].X / gas2['O2'].X) / (4.0*ct.faraday))

def Solver(f, xstart, C=0.0):
    func = lambda x:f(x)-C
    return fsolve(func, 0.0)

# def Solver(f, xstart, C=0.0):
#     f0 = f(xstart) - C
#     x0 = xstart
#     dx = 1.0e-6
#     n = 0
#     while n < 200:
#         ff = f(x0 + dx) - C
#         dfdx = (ff - f0)/dx
#         step = - f0/dfdx

#         # avoid taking steps too large
#         if abs(step) > 0.1:
#             step = 0.1*step/abs(step)

#         x0 += step
#         emax = 0.00001  # 0.01 mV tolerance
#         if abs(f0) < emax and n > 8:
#             return x0
#         f0 = f(x0) - C
#         n += 1
#     raise Exception('no root!')
   
if __name__=="__main__":
    # test values
    print('[ 1.30661931e+03  2.45000000e-01 -5.23098811e-01  3.26654827e-01  4.28440650e-02]')
    print('[ 5.82664760e+02  2.45000000e-01 -5.66693531e-01  1.45666190e-01  1.90306283e-01]')
    # output must be the same
    print(TSOFC(T=1073.0).getPolarizationCurve()[-1])
    print(TSOFC(T=1000.0).getPolarizationCurve()[-1])

    # import matplotlib
    # matplotlib.use('Qt5Agg')
    # import matplotlib.pyplot as plt   
    # fig = plt.figure()
    # ax = fig.gca()
    # plot_data = TSOFC(T=1073.0).getPolarizationCurve()
    # ax.plot(plot_data[:,0],plot_data[:,-1],"-")
    # ax.set_xlabel(r"I, [A]")
    # ax.set_ylabel(r"U, [В]")
    # # plt.savefig("tests/tredox_polcurve.pdf")
    # plt.show()
