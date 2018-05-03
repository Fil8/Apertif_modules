#!/usr/bin/env python

# Import modules
import numpy as np
import pyrap.tables as tables
import sys



############ START USER SETTINGS ############

# Input file (same channelisation required, different array layouts allowed)
MS=[
    '/scratch/paolo/all_ic5264/msdir/1491463063.ms',
    '/scratch/paolo/all_ic5264/msdir/1491480644.ms',
   ]

# Select field for which to calculate the theoretical noise
selectFieldName='IC5264'         # set to None to select all fields

# Plot?
makePlot=True

# Single antenna specs (the number of antennas in the array is derived from the input .MS files)
k = 1380.6                   # Jy m^2 / K
Tsys=22                      # K
eff=1.                       # antenna efficiency
diam=13.5                    # m
Aant=np.pi*(diam/2)**2       # collecting area of 1 antenna
SEFD=2*k*Tsys/eff/Aant       # Jy

############ END USER SETTINGS ############



# Define functions
def getData(ms,k_value,Tsys_value,eff_value,Aant_value):
    print ''
    print '--- Working on file {0:s} ---'.format(ms)
    t=tables.table(ms)
    fieldIDs=t.getcol('FIELD_ID')
    ant1=t.getcol('ANTENNA1')
    ant2=t.getcol('ANTENNA2')
    fieldNames=tables.table(ms+'/FIELD').getcol('NAME')
    spw=tables.table(ms+'/SPECTRAL_WINDOW')
    channelWidths=spw.getcol('CHAN_WIDTH')
    channelFreqs=spw.getcol('CHAN_FREQ')

    if selectFieldName:
        try:
            selectFieldID=fieldNames.index(selectFieldName)
        except ValueError:
            print ' CATASTROPHE!'
            print ' Cannot find the field you want to process, {0:s}'.format(selectFieldName)
            print ' Available fields are',fieldNames
            print ' Aborting ...'
            exit()
        print 'Successfully selected Field with name {0:s} (Field ID = {1:d})'.format(selectFieldName,selectFieldID)
        selection=fieldIDs==selectFieldID
    else:
        print ' Will process all available fields:',fieldNames
        selection=fieldIDs>=fieldIDs.min()

    autoCorr=ant1==ant2
    if autoCorr.sum(): print 'Successfully selected crosscorrelations only'
    else: print 'Found crosscorrelations only'
    selection*=ant1!=ant2
    nrAnt=np.unique(np.concatenate((ant1,ant2))).shape[0]
    nrBaseline=nrAnt*(nrAnt-1)/2
    print 'Number of antennas  = {0:d}'.format(nrAnt)
    print 'Number of baselines = {0:d}'.format(nrBaseline)
    print 'Frequency coverage  = {0:.5e} Hz - {1:.5e} Hz'.format(channelFreqs.min(),channelFreqs.max())
    if np.unique(channelWidths).shape[0]==1: print 'Channel width = {0:.5e} Hz'.format(np.unique(channelWidths)[0])
    else: print 'The channel width takes the following values:',channelWidths,'Hz'

    print 'Loading flags and intervals ...'
    flag=t.getcol('FLAG')[selection]          # flagged data have flag = True
    interval=t.getcol('INTERVAL')[selection]
    t.close()

    print 'The *flag* array has shape (Nr_integrations, Nr_channels, Nr_polarisations) =',flag.shape
    print 'The *interval* array has shape (Nr_integrations) =',interval.shape
    print 'The *channel* width array has shape (-, Nr_channels) =',channelWidths.shape

    rms=np.sqrt(2)*k_value*Tsys_value/eff_value/Aant_value/np.sqrt(channelWidths*interval.sum()*flag.shape[2])
    if len(rms.shape)==2 and rms.shape[0]==1: rms=rms[0]

    print 'Total Integration on selected field(s) = {0:.2f} h ({1:d} polarisations)'.format(interval.sum()/nrBaseline/3600,flag.shape[2])
    print 'The Stokes I theoretical rms ignoring flags is in the range:    *** ({0:.3e} - {1:.3e}) Jy ***'.format(np.nanmin(rms),np.nanmax(rms))

    return flag,interval,channelWidths,channelFreqs,rms


# Print assumptions
print ''
print '--- Assumptions ---'
print '  Tsys/efficiency      = {0:.1f} K (frequency independent)'.format(Tsys/eff)
print '  Dish diameter        = {0:.1f} m'.format(diam)
print '    and therefore SEFD = {0:.1f} Jy'.format(SEFD)


# Read MS files to get the flags and calculate single-file rms

# Start with first file ...
flag0,interval0,channelWidths0,channelFreqs0,rms0=getData(MS[0],k,Tsys,eff,Aant)
rmsAll=[rms0]

# ... and do the same for all other MS's appending to the flag array, checking that the channelisation is the same
for ii in range(1,len(MS)):
    flagi,intervali,channelWidthsi,channelFreqsi,rmsi=getData(MS[ii],k,Tsys,eff,Aant)

    if (channelWidths0!=channelWidthsi).sum() or (channelFreqs0!=channelFreqsi).sum():
        print ' CATASTROPHE!'
        print ' The {0:d}-th input .MS file {1:s} has different channelization than the 0-th input .MS file {2:s}'.format(ii,MS[ii],MS[0])
        print ' Cannot combine files to estimate their joint theoretical noise'
        print ' Aborting ...'
        exit()
    else:
        flag0=np.concatenate((flag0,flagi),axis=0)
        interval0=np.concatenate((interval0,intervali),axis=0)
        rmsAll.append(rmsi)


# Concatenate files
print ''
print '--- All input tables concatenated ---'
print 'The concatenated *flag* array has shape (Nr_integrations, Nr_channels, Nr_polarisations) =',flag0.shape
print 'The concatenated *interval* array has shape (Nr_integrations) =',interval0.shape
print 'The concatenated *channel* width array has shape (-, Nr_channels) =',channelWidths0.shape

print ''
print '--- Result ---'
print 'Calculating rms ...'


# Reshape arrays
interval0.resize((interval0.shape[0],1,1))
channelWidths0.resize((channelWidths0.shape[1]))
channelFreqs0.resize((channelFreqs0.shape[1]))


# Calculate theoretical rms
rmsAll=np.array(rmsAll)
rmsAll=1./np.sqrt( (1./rmsAll**2).sum(axis=0) )
unflaggedIntegration=(interval0*(1-flag0.astype(int))).sum(axis=(0,2)) # total integration per channel adding up all UNFLAGGED integrations and polarisations (sec)
unflaggedIntegration[unflaggedIntegration==0]=np.nan
rmsUnflagged=np.sqrt(2)*k*Tsys/eff/Aant/np.sqrt(channelWidths0*unflaggedIntegration)

print 'The Stokes I theoretical rms ignoring flags is in the range:    *** ({0:.3e} - {1:.3e}) Jy ***'.format(np.nanmin(rmsAll),np.nanmax(rmsAll))
print 'The Stokes I theoretical rms applying flags is in the range:    *** ({0:.3e} - {1:.3e}) Jy ***'.format(np.nanmin(rmsUnflagged),np.nanmax(rmsUnflagged))


# Plot
if makePlot==True:
    print ''
    print '--- Plot ---'
    print 'Busy plotting ...'
    import matplotlib.pyplot as plt
    plt.subplot(111)
    plt.plot(channelFreqs0,rmsAll*1e+3,'k:')
    plt.plot(channelFreqs0,rmsUnflagged*1e+3,'r-')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('rms (mJy)')
    plt.legend(('all data','unflagged data'),loc='lower right')
    plt.savefig('rms.png')
    print 'rms.png saved in working directory'
