import matplotlib.pyplot as plt
import numpy as np
import copy
import pickle
import matplotlib

# output figure filename
filename = "stacked_plot"
molecule_name = "noradrenalin"

# Directories to look for data
datadirs = ['noradrenalin_zeropadded']
linewidths = [1.2]


fig1 = None
ax11 = None
ax12 = None
ax21 = None
ax22 = None
leg_real = None
leg_imag = None
leg_field = None

def mirror_and_zeropad(arr, zeropad):
    arr_inv = copy.copy(arr[1:][::-1])
    arr_inv = [np.conj(xv) for xv in arr_inv]
    arr = arr_inv + arr
    arr = zeropad+arr+zeropad
    return arr


for i, (lw, directory) in enumerate(zip(linewidths, datadirs)):
    field_w = np.load('{}/{}_field_w.npy'.format(directory, molecule_name))
    freq = np.load('{}/{}_freq.npy'.format(directory, molecule_name))
    f = open('{}/prt_{}_function.pickle'.format(directory, molecule_name), 'rb')
    prop = pickle.load(f)
    properties = {}
  
    # load pulse settings
    pulse_settings = prop['pulse_settings'] if 'pulse_settings' in prop else {}

    # Use zero padded numbers if available
    if 'properties_zeropad' in prop:
        properties = prop['properties_zeropad']
    else:
        properties = prop['properties'] if 'properties' in prop else prop
    
    xdata = []
    ydata = []
    zdata = []
    for key, item in properties.items():
        _freq = key[2]
        key = (str(key[0]), str(key[1]))

        # Debug Veloxchem to verify the new factor -1
        item *= -1.
        if key == ('x', 'x'):
            xdata.append((_freq, item))
        if key == ('y','y'):
            ydata.append((_freq, item))
        if key == ('z','z'):
            zdata.append((_freq, item))

    
    xdata.sort(key=lambda pair: pair[0])
    ydata.sort(key=lambda pair: pair[0])
    zdata.sort(key=lambda pair: pair[0])
   
    # Extract frequencies
    xfreqs = [data[0] for data in xdata]
    yfreqs = [data[0] for data in ydata]
    zfreqs = [data[0] for data in zdata]
    x_alphas = [data[1] for data in xdata]
    y_alphas = [data[1] for data in ydata]
    z_alphas = [data[1] for data in zdata]

    # Add mirror
    xfrec_inv = copy.copy(xfreqs[1:][::-1])
    xfrec_inv = [xf*-1. for xf in xfrec_inv]
    xfreqs = xfrec_inv + xfreqs
    
    # Zero padding
    zeropad = [0.0] * 5000

    # Add mirror and zeropadding
    x_alphas = mirror_and_zeropad(x_alphas, zeropad)
    y_alphas = mirror_and_zeropad(y_alphas, zeropad)
    z_alphas = mirror_and_zeropad(z_alphas, zeropad)

    # Add mirror
    field_w = field_w.tolist()
    field_w_inv = copy.copy(field_w[1:][::-1])
    field_w_inv = [np.conj(xv) for xv in field_w_inv]
    field_w = field_w_inv + field_w
    
    field_w = zeropad+field_w+zeropad
    
    dxfreq = xfreqs[1] - xfreqs[0]
    dfreq = dxfreq * (len(x_alphas)-1)/2.
    xfreqs = np.linspace(-dfreq, dfreq, len(field_w), endpoint=True)
    
    field_w = np.array(field_w)
    x_alphas = np.array(x_alphas)
    y_alphas = np.array(y_alphas)
    z_alphas = np.array(z_alphas)
    xfreqs = np.array(xfreqs)

    # Next step is to try generate the times based on the frequencies read in
    w = xfreqs
    dw=w[1] - w[0]
    N = len(xfreqs)
    dt = 2.*np.pi / (N * dw)
    t_end = ((N-1)/2) * dt
    t = np.linspace(-t_end, t_end, len(xfreqs), endpoint=True)

    tshift = np.fft.ifftshift(t) # 0, positive ascending, negative ascending
    wshift = np.fft.ifftshift(w) # 0, positive ascending, negative ascending
    

    # start by shift alpha to the same order as Fw, i.e., 0, positive ascending, negative ascending and then perform DFT
    # Multiply by phase factor to get continuous FT
    alphas = 1./3. * (x_alphas+y_alphas+z_alphas) # To get isotropic
    dipmom = []
    dipmom = np.fft.fft(np.fft.ifftshift(alphas*field_w)) 
    dipmom *= dw*np.exp(-1.j*wshift[0]*tshift)
    dipmom_plot = np.fft.fftshift(dipmom)

    # Calculate the time series of the field
    field_t = []
    field_t = np.fft.fft(np.fft.ifftshift(field_w))
    field_t *= dw*np.exp(-1.j*wshift[0]*tshift)
    field_t_plot = np.fft.fftshift(field_t)

    # -------------------------------------------
    # 	Time series plot
    # -------------------------------------------
    # Initialize figure and axis if first loop
    if fig1 is None:
       fig1 = plt.figure(1, figsize=(10,5))

    if ax11 is None: 
       (ax21, ax11) = fig1.subplots(2,1)

    if ax12 is None: ax12 = ax11.twinx()

    ax12.plot(t,np.real(dipmom_plot)*1.e4,color='r', linewidth=lw)
    ax12.plot(t,np.imag(dipmom_plot)*1.e4,color='r',linestyle='--', linewidth=lw)

    # Plotting field
    ax11.plot(t,np.real(field_t_plot)*1e5,color='k', linewidth=1.2)

    # Adding labels
    ax11.text(0.01, 1.01, r'$\times 10^{-5}$', transform=ax11.transAxes)
    ax11.text(0.94, 1.01, r'$\times 10^{-4}$', transform=ax11.transAxes)
    ax11.set_xlabel('Time (a.u.)')
    ax11.set_ylabel(r'$F$(t) (a.u.)')
    ax12.set_ylabel(r'$\mu$(t) (a.u.)',color='r')

    # Calculate the plot limits
    dipmax = np.max((np.abs(np.real(dipmom_plot)))) * 1.e4
    fieldmax = np.max((np.abs(np.real(field_t_plot)))) *1.e5 

    # Set the plot limits
    ax11.set_ylim([-fieldmax,fieldmax])
    ax12.set_ylim([-dipmax,dipmax])
    ax11.set_xlim([-max(t), max(t)])
    
    
    # -------------------------------------------
    # 	Frequency plot
    # -------------------------------------------
    if ax22 is None: ax22 = ax21.twinx()

    leg_real, = ax22.plot(xfreqs,np.real(alphas), color='r', linewidth=lw, label='real')
    leg_imag, = ax22.plot(xfreqs,np.imag(alphas), color='b', linewidth=lw, label='imag')

    # Plotting field
    leg_field, = ax21.plot(xfreqs, np.real(field_w)*1.e5, color='k', linestyle='-', lw=0.4, label='field')
    ax21.plot(xfreqs, np.abs(field_w)*1.e5, color='k', linestyle=':', lw=0.4)

    # Add label
    ax21.text(0.01, 1.01, r'$\times 10^{-5}$', transform=ax21.transAxes)

    max_field = np.max(np.real(field_w))
    max_oscillator_strengths = np.max(prop['rsp_results']['oscillator_strengths'])

    # Plot sticks from eigenvalues
    prop['rsp_results']['oscillator_strengths'] = [x/max_oscillator_strengths * max_field for x in prop['rsp_results']['oscillator_strengths']]
    for fre, osc in zip(prop['rsp_results']['eigenvalues'], prop['rsp_results']['oscillator_strengths']):
        ax21.plot([fre, fre],[0, osc], color='k', lw=1)

    # Set the plot limits
    x_alphamax = max(x_alphas) * 1.7
    ax22.set_ylim([-x_alphamax, x_alphamax])
    ax21.set_xlim([-1, 1])

# ---------------------------------------------
#   Add legend 
# ---------------------------------------------
ax21.legend([leg_real, leg_imag, leg_field], ['real', 'imag','field'], frameon=False)

# ---------------------------------------------
#   Save the plot
# ---------------------------------------------
ax21.set_ylabel(r'F($\omega$)  (a.u.)')
ax22.set_ylabel(r'$\overline{\alpha}$($\omega$)  (a.u.)')
ax21.set_xlabel('Frequencies (a.u.)')
ax11.set_xlabel('Time (a.u.)')

plt.subplots_adjust(hspace=0.34, top=0.96)

plt.show()
fig1.savefig('{}.pdf'.format(filename),dpi=300)
   
plt.close(fig1)

