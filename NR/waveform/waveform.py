import bilby

from numpy import sqrt, zeros, pi, cos, sin, linspace
import copy
from astropy.cosmology import Planck15, z_at_value
import astropy.units as u
from matplotlib import pyplot as plt
z=1.0
d_L = Planck15.luminosity_distance(z)
print(d_L)
d_L = d_L.value

params_wf_gen = dict(
            mass_1              = 1.4*(1+z), 
            mass_2              = 1.4*(1+z),
            luminosity_distance = d_L, 
            theta_jn            = 0.0, 
            phase               = 0.0, 
            geocent_time        = 0.0, 
            a_1                 = 0.0, 
            a_2                 = 0.0, 
            tilt_1              = 0.0, 
            tilt_2              = 0.0,
            phi_12              = 0.0, 
            phi_jl              = 0.0)

waveform_arguments = dict(
            waveform_approximant='IMRPhenomD',
            minimum_frequency=20.0)
waveform_generator = bilby.gw.WaveformGenerator(
            sampling_frequency=2048,
            duration=128,
            frequency_domain_source_model = bilby.gw.source.lal_binary_black_hole,
            parameters = params_wf_gen,
            parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
            waveform_arguments = waveform_arguments)
f_array = waveform_generator.frequency_array
h_plus = waveform_generator.frequency_domain_strain()['plus']
h_cross = waveform_generator.frequency_domain_strain()['cross']
# print(f_array)
# print(h_plus)
# print(h_cross)

theta = 0.0
phi = 0.0
psi = 1.0
F_plus  = sqrt(3)/2 * ( 1/2*(1+cos(theta)**2)*cos(2*phi)*cos(2*psi) - cos(theta)*sin(2*phi)*sin(2*psi) )
F_cross = sqrt(3)/2 * ( 1/2*(1+cos(theta)**2)*cos(2*phi)*sin(2*psi) + cos(theta)*sin(2*phi)*cos(2*psi) )
# print(F_plus)
# print(F_cross)

hf = F_plus*h_plus + F_cross*h_cross
# print(hf)

def loadData(filename):
    inFile=open(filename,'r')
    y=[]
    for line in inFile:
        y.append(float(line))
    return y

plt.figure()
plt.semilogx(f_array, abs(hf),'r')
plt.title('abs(hf)')
plt.xlim(20,1000)
#plt.savefig('abshf.png')
y=loadData('H.txt')
length=len(y)
x=linspace(20,1000,length)
plt.semilogx(x,y,'b')
plt.show()

plt.figure()
plt.semilogx(f_array, hf.real,'r')
plt.title('hf_real')
plt.xlim(20,1000)
#plt.savefig('hf_real.png')
y=loadData('hf_real.txt')
length=len(y)
x=linspace(20,1000,length)
plt.semilogx(x,y,'b')
plt.show()
    #plt.savefig('hf_real.png')


#plt.figure()
#plt.semilogx(f_array, hf.imag)
#plt.title('hf_imag')
#plt.xlim(20,1000)
#plt.savefig('hf_imag.png')

#plt.figure()
#plt.semilogx(f_array, h_plus.imag)
#plt.title('h_plus_imag')
#plt.xlim(20,1000)
#plt.savefig('h_plus_imag.png')

#plt.figure()
#plt.semilogx(f_array, h_cross.imag)
#plt.title('h_cross_imag')
#plt.xlim(20,1000)
#plt.savefig('h_cross_imag.png')














# ################################################################################
# # define the PV source model
# def waveform_PV(
#         frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
#         phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, PV_A, **kwargs):
#     waveform_GR = bilby.gw.source.lal_binary_black_hole(
#         frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2, 
#         luminosity_distance=luminosity_distance, a_1=a_1, tilt_1=tilt_1,
#         phi_12=phi_12, a_2=a_2, tilt_2=tilt_2, phi_jl=phi_jl, theta_jn=theta_jn,
#         phase=phase, **kwargs)
#     h_plus  = waveform_GR['plus']
#     h_cross = waveform_GR['cross']
    
#     frequency_kwargs = dict(minimum_frequency=20.0, 
#                             maximum_frequency=frequency_array[-1])
#     frequency_kwargs.update(kwargs)
#     minimum_frequency = frequency_kwargs['minimum_frequency']
#     maximum_frequency = frequency_kwargs['maximum_frequency']
#     frequency_bounds = ((frequency_array >= minimum_frequency) * 
#                         (frequency_array <= maximum_frequency))

#     f_array = frequency_array * frequency_bounds
#     dphi1 = PV_A*(pi*f_array)**2
#     h_plus_PV = h_plus + h_cross*dphi1
#     h_cross_PV = h_cross - h_plus*dphi1

#     return dict(plus=h_plus_PV, cross=h_cross_PV)
# ################################################################################

# def fisher_matrix(detector_name, fiducial_dict, step_dict):
#     ifo = bilby.gw.detector.get_empty_interferometer(detector_name)
#     if type(ifo) == bilby.gw.detector.networks.TriangularInterferometer:
#         f_matrix_ls = []
#         SNR2_ls = []
#         for item in ifo:
#             f_matrix_i, SNR2_i = fisher_matrix_for_single_ifo(item, fiducial_dict, step_dict)
#             f_matrix_ls.append(f_matrix_i)
#             SNR2_ls.append(SNR2_i)
#         f_matrix = sum(f_matrix_ls)
#         SNR2 = sum(SNR2_ls)
#     elif type(ifo) == bilby.gw.detector.interferometer.Interferometer:
#         f_matrix, SNR2 = fisher_matrix_for_single_ifo(ifo, fiducial_dict, step_dict)
#     else:
#         raise TypeError('{} is not an expected object'.format(type(ifo)))

#     return f_matrix, SNR2        


# def fisher_matrix_for_single_ifo(ifo, fiducial_dict, step_dict):
#     '''
#     compute Fisher matrix for a interferometer
#     note: total mass and d_L are directly passed to bilby.WaveformGenerator.
#           considering redshift mass and scaling d_L in main program.

#     Parameters:
#     --------
#     ifo: bilby.Interferometer object
#         returned by bilby.gw.detector.get_empty_interferometer
#     fiducial_dict: dict
#         fiducial parameters dict. 
#         total_mass(solar mass), eta, d_L(Mpc), 
#         ra(radians), dec(radians), iota(radians), t_c(geocentric time), 
#         ph_c(radians), pol_angle(radians),
#         PV_A(s^2)
#     step_dict: dict
#         step of each parameter

#     Retures:
#     --------
#     Fisher matrix: array

#     '''
#     # setting
#     # sampling_frequency = 16384
#     # sampling_frequency = 1024
#     # sampling_frequency = 2048
#     sampling_frequency = 4096
#     # sampling_frequency = 8192
#     duration = 4
#     # waveform_approximant='TaylorF2'
#     waveform_approximant='IMRPhenomPv2'

#     # initializing a Interferometer object 
#     # getting frequency_array and PSD_array
#     # ifo = bilby.gw.detector.get_empty_interferometer(detector_name)
#     ifo.set_strain_data_from_zero_noise(
#         sampling_frequency=sampling_frequency, 
#         duration=duration)
#     f_array = ifo.frequency_array
#     minimum_frequency = ifo.minimum_frequency
#     # f_mask = ifo.frequency_mask
#     PSD = ifo.power_spectral_density
#     # PSD_array = ifo.power_spectral_density_array

#     # hf the strain data of the detector
#     def hf(parameters):
#         M   = parameters['total_mass']
#         eta = parameters['eta']
#         params_wf_gen = dict(
#             mass_1              = M*(1 + sqrt(1-4*eta))/2, 
#             mass_2              = M*(1 - sqrt(1-4*eta))/2,
#             ra                  = parameters['ra'], 
#             dec                 = parameters['dec'], 
#             luminosity_distance = parameters['d_L'], 
#             theta_jn            = parameters['iota'], 
#             psi                 = parameters['pol_angle'],
#             phase               = parameters['ph_c'], 
#             geocent_time        = parameters['t_c'], 
#             a_1                 = 0.0, 
#             a_2                 = 0.0, 
#             tilt_1              = 0.0, 
#             tilt_2              = 0.0,
#             phi_12              = 0.0, 
#             phi_jl              = 0.0,
#             PV_A                = parameters['PV_A'])
#         waveform_arguments = dict(
#             waveform_approximant=waveform_approximant,
#             minimum_frequency=minimum_frequency)
#         waveform_generator = bilby.gw.WaveformGenerator(
#             sampling_frequency=sampling_frequency,
#             duration=duration,
#             frequency_domain_source_model = waveform_PV,
#             parameters = params_wf_gen,
#             parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
#             waveform_arguments = waveform_arguments)
#         hf = ifo.get_detector_response(
#             waveform_generator.frequency_domain_strain(), 
#             params_wf_gen)
#         return hf
    
#     # generate Fisher matrix
#     params_d_total_mass = copy.deepcopy(fiducial_dict)
#     params_d_eta        = copy.deepcopy(fiducial_dict)
#     params_d_d_L        = copy.deepcopy(fiducial_dict)
#     params_d_ra         = copy.deepcopy(fiducial_dict)
#     params_d_dec        = copy.deepcopy(fiducial_dict)
#     params_d_iota       = copy.deepcopy(fiducial_dict)
#     params_d_t_c        = copy.deepcopy(fiducial_dict)
#     params_d_ph_c       = copy.deepcopy(fiducial_dict)
#     params_d_pol_angle  = copy.deepcopy(fiducial_dict)
#     params_d_PV_A       = copy.deepcopy(fiducial_dict)
#     params_d_total_mass['total_mass'] += step_dict['total_mass']
#     params_d_eta['eta']               += step_dict['eta']
#     params_d_d_L['d_L']               += step_dict['d_L']
#     params_d_ra['ra']                 += step_dict['ra']
#     params_d_dec['dec']               += step_dict['dec']
#     params_d_iota['iota']             += step_dict['iota']
#     params_d_t_c['t_c']               += step_dict['t_c']
#     params_d_ph_c['ph_c']             += step_dict['ph_c']
#     params_d_pol_angle['pol_angle']   += step_dict['pol_angle']
#     params_d_PV_A['PV_A']             += step_dict['PV_A']
#     hf_fidu               = hf(fiducial_dict       )
#     hf_d_total_mass       = hf(params_d_total_mass ) 
#     hf_d_eta              = hf(params_d_eta        )
#     hf_d_d_L              = hf(params_d_d_L        )
#     hf_d_ra               = hf(params_d_ra         )
#     hf_d_dec              = hf(params_d_dec        )
#     hf_d_iota             = hf(params_d_iota       )
#     hf_d_t_c              = hf(params_d_t_c        )
#     hf_d_ph_c             = hf(params_d_ph_c       )
#     hf_d_pol_angle        = hf(params_d_pol_angle  )
#     hf_d_PV_A             = hf(params_d_PV_A       )

#     dh_lst = [ (hf_d_total_mass - hf_fidu)/step_dict['total_mass'],
#                (hf_d_eta        - hf_fidu)/step_dict['eta'],
#                (hf_d_d_L        - hf_fidu)/step_dict['d_L'],
#                (hf_d_ra         - hf_fidu)/step_dict['ra'],
#                (hf_d_dec        - hf_fidu)/step_dict['dec'],
#                (hf_d_iota       - hf_fidu)/step_dict['iota'],
#                (hf_d_t_c        - hf_fidu)/step_dict['t_c'],
#                (hf_d_ph_c       - hf_fidu)/step_dict['ph_c'],
#                (hf_d_pol_angle  - hf_fidu)/step_dict['pol_angle'],
#                (hf_d_PV_A       - hf_fidu)/step_dict['PV_A'] ]
#     num_para = len(dh_lst)
#     F_matrix = zeros((num_para, num_para))
#     for i in range(0, num_para):
#         for j in range(0, i+1):
#             m_e = bilby.gw.utils.inner_product(dh_lst[i], dh_lst[j], frequency=f_array, PSD=PSD)
#             F_matrix[i][j] = m_e
#             F_matrix[j][i] = m_e

#     # snr
#     SNR2 = bilby.gw.utils.inner_product(hf_fidu, hf_fidu, frequency=f_array, PSD=PSD)
    
#     return F_matrix, SNR2
