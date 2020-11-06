import bilby
import logging
from numpy import sqrt, zeros, pi
from numpy.linalg import inv
import copy
from astropy.cosmology import Planck15, z_at_value
import astropy.units as u



def fisher_matrix(detector_name, fiducial_dict, step_dict):
    '''
    allowed detector_name: 'L1', 'H1', 'V1', 'K1', 'CE', 'ET', 'GEO600'
    '''
    ifo = bilby.gw.detector.get_empty_interferometer(detector_name)
    if type(ifo) == bilby.gw.detector.networks.TriangularInterferometer:
        f_matrix_ls = []
        SNR2_ls = []
        for item in ifo:
            f_matrix_i, SNR2_i = fisher_matrix_for_single_ifo(item, fiducial_dict, step_dict)
            f_matrix_ls.append(f_matrix_i)
            SNR2_ls.append(SNR2_i)
        f_matrix = sum(f_matrix_ls)
        SNR2 = sum(SNR2_ls)
    elif type(ifo) == bilby.gw.detector.interferometer.Interferometer:
        f_matrix, SNR2 = fisher_matrix_for_single_ifo(ifo, fiducial_dict, step_dict)
    else:
        raise TypeError('{} is not an expected object'.format(type(ifo)))

    return f_matrix, SNR2        


def fisher_matrix_for_single_ifo(ifo, fiducial_dict, step_dict):
    '''
    compute Fisher matrix for a interferometer
    note: total mass and d_L are directly passed to bilby.WaveformGenerator.
          considering redshift mass and scaling d_L in main program.

    Parameters:
    --------
    ifo: bilby.Interferometer object
        returned by bilby.gw.detector.get_empty_interferometer
    fiducial_dict: dict
        fiducial parameters dict. 
        total_mass(solar mass), eta, d_L(Mpc), 
        ra(radians), dec(radians), iota(radians), t_c(geocentric time), 
        ph_c(radians), pol_angle(radians),
        PV_A(s^2)
    step_dict: dict
        step of each parameter

    Retures:
    --------
    Fisher matrix: array

    '''
    # setting
    # sampling_frequency = 16384
    # sampling_frequency = 1024
    # sampling_frequency = 2048
    sampling_frequency = 4096
    # sampling_frequency = 8192
    duration = 4
    # waveform_approximant='TaylorF2'
    waveform_approximant='IMRPhenomPv2'

    # initializing a Interferometer object 
    # getting frequency_array and PSD_array
    # ifo = bilby.gw.detector.get_empty_interferometer(detector_name)
    ifo.set_strain_data_from_zero_noise(
        sampling_frequency=sampling_frequency, 
        duration=duration)
    f_array = ifo.frequency_array
    minimum_frequency = ifo.minimum_frequency
    # f_mask = ifo.frequency_mask
    PSD = ifo.power_spectral_density
    # PSD_array = ifo.power_spectral_density_array

    # hf the strain data of the detector
    def hf(parameters):
        M   = parameters['total_mass']
        eta = parameters['eta']
        params_wf_gen = dict(
            mass_1              = M*(1 + sqrt(1-4*eta))/2, 
            mass_2              = M*(1 - sqrt(1-4*eta))/2,
            ra                  = parameters['ra'], 
            dec                 = parameters['dec'], 
            luminosity_distance = parameters['d_L'], 
            theta_jn            = parameters['iota'], 
            psi                 = parameters['pol_angle'],
            phase               = parameters['ph_c'], 
            geocent_time        = parameters['t_c'], 
            a_1                 = 0.0, 
            a_2                 = 0.0, 
            tilt_1              = 0.0, 
            tilt_2              = 0.0,
            phi_12              = 0.0, 
            phi_jl              = 0.0)
        waveform_arguments = dict(
            waveform_approximant=waveform_approximant,
            minimum_frequency=minimum_frequency)
        waveform_generator = bilby.gw.WaveformGenerator(
            sampling_frequency=sampling_frequency,
            duration=duration,
            frequency_domain_source_model = bilby.gw.source.lal_binary_black_hole,
            parameters = params_wf_gen,
            parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
            waveform_arguments = waveform_arguments)
        hf = ifo.get_detector_response(
            waveform_generator.frequency_domain_strain(), 
            params_wf_gen)
        return hf
    
    # generate Fisher matrix
    params_d_total_mass = copy.deepcopy(fiducial_dict)
    params_d_eta        = copy.deepcopy(fiducial_dict)
    params_d_d_L        = copy.deepcopy(fiducial_dict)
    # params_d_ra         = copy.deepcopy(fiducial_dict)
    # params_d_dec        = copy.deepcopy(fiducial_dict)
    params_d_iota       = copy.deepcopy(fiducial_dict)
    params_d_t_c        = copy.deepcopy(fiducial_dict)
    params_d_ph_c       = copy.deepcopy(fiducial_dict)
    params_d_pol_angle  = copy.deepcopy(fiducial_dict)
    params_d_total_mass['total_mass'] += step_dict['total_mass']
    params_d_eta['eta']               += step_dict['eta']
    params_d_d_L['d_L']               += step_dict['d_L']
    # params_d_ra['ra']                 += step_dict['ra']
    # params_d_dec['dec']               += step_dict['dec']
    params_d_iota['iota']             += step_dict['iota']
    params_d_t_c['t_c']               += step_dict['t_c']
    params_d_ph_c['ph_c']             += step_dict['ph_c']
    params_d_pol_angle['pol_angle']   += step_dict['pol_angle']
    hf_fidu               = hf(fiducial_dict       )
    hf_d_total_mass       = hf(params_d_total_mass ) 
    hf_d_eta              = hf(params_d_eta        )
    hf_d_d_L              = hf(params_d_d_L        )
    # hf_d_ra               = hf(params_d_ra         )
    # hf_d_dec              = hf(params_d_dec        )
    hf_d_iota             = hf(params_d_iota       )
    hf_d_t_c              = hf(params_d_t_c        )
    hf_d_ph_c             = hf(params_d_ph_c       )
    hf_d_pol_angle        = hf(params_d_pol_angle  )

    dh_lst = [ (hf_d_total_mass - hf_fidu)/step_dict['total_mass'],
               (hf_d_eta        - hf_fidu)/step_dict['eta'],
               (hf_d_t_c        - hf_fidu)/step_dict['t_c'],
               (hf_d_ph_c       - hf_fidu)/step_dict['ph_c'],
               (hf_d_iota       - hf_fidu)/step_dict['iota'],
               (hf_d_pol_angle  - hf_fidu)/step_dict['pol_angle'],
               (hf_d_d_L        - hf_fidu)/step_dict['d_L']
            #    (hf_d_ra         - hf_fidu)/step_dict['ra'],
            #    (hf_d_dec        - hf_fidu)/step_dict['dec'],
            ]
    num_para = len(dh_lst)
    F_matrix = zeros((num_para, num_para))
    for i in range(0, num_para):
        for j in range(0, i+1):
            m_e = bilby.gw.utils.inner_product(dh_lst[i], dh_lst[j], frequency=f_array, PSD=PSD)
            F_matrix[i][j] = m_e
            F_matrix[j][i] = m_e

    # snr
    SNR2 = bilby.gw.utils.inner_product(hf_fidu, hf_fidu, frequency=f_array, PSD=PSD)
    
    return F_matrix, SNR2


################################################################################

total_mass= 3.0
d_L       = 475.0
eta       = 0.24
ra        = pi
dec       = 0.
iota      = pi/2
t_c       = 1126259462
ph_c      = pi/4
pol_angle = pi/3
fiducial_dict = dict(
    total_mass = total_mass,
    d_L        = d_L,
    eta        = eta,
    ra         = ra , 
    dec        = dec , 
    iota       = iota , 
    t_c        = t_c , 
    ph_c       = ph_c , 
    pol_angle  = pol_angle)
step_dict = dict(    
    total_mass = 1e-8, 
    eta        = 1e-8, 
    d_L        = 1e-5, 
    ra         = 1e-6, 
    dec        = 1e-6, 
    iota       = 1e-4, 
    t_c        = 1e-6, 
    ph_c       = 1e-3, 
    pol_angle  = 1e-3)



fm_CE, SNR2_CE = fisher_matrix('CE', fiducial_dict, step_dict)
fm_ET, SNR2_ET = fisher_matrix('ET', fiducial_dict, step_dict)

fm = fm_ET
SNR2 = SNR2_ET

print(sqrt(SNR2))
print(fm)
cov_m = inv(fm)
print(cov_m)
print('*********************')
print(sqrt(cov_m[0,0]))
print(sqrt(cov_m[1,1]))
print(sqrt(cov_m[2,2]))
print(sqrt(cov_m[3,3]))
print(sqrt(cov_m[4,4]))
print(sqrt(cov_m[5,5]))
print(sqrt(cov_m[6,6])/d_L)
# print(sqrt(cov_m[7,7]))
# print(sqrt(cov_m[8,8]))
print('*********************')
print(fm[0,0])
print(fm[1,1])
print(fm[2,2])
print(fm[3,3])
print(fm[4,4])
print(fm[5,5])
print(fm[6,6])
# print(fm[7,7])
# print(fm[8,8])




