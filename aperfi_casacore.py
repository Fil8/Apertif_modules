#!/usr/bin/env python
__author__ = "Filippo Maccagni, Tom Osterloo"
__copyright__ = "ASTRON"
__email__ = "maccagni@astron.nl"

# Import modules
import numpy as np
import pyrap.tables as tables
import numpy.ma as ma
import aplpy
import astropy.io.fits as pyfits

############ START USER SETTINGS ############

# Input file (same channelisation required, different array layouts allowed)
MS=['/home/maccagni/Projects/ApertiF/RFI/data/WSRTA110311_B000.MS']
   
############ END USER SETTINGS ############
def write_freq_base(rms,channelWidths,channelFreqs,fitsfile) :
	'''
	Writes an image.fits of visibilities ordered by frequency, baseline.
	Baselines are ordered by their length.
	'''
	#reverse frequencies if going from high-to-low         
	
	#set fits file
	hdu = pyfits.PrimaryHDU(rms)
	hdulist = pyfits.HDUList([hdu])
	header = hdulist[0].header
	print channelWidths[0,0]
	print rms.shape
	print channelFreqs[0,0]
	#write header keywords               
	header['CRPIX1'] = 1
	header['CDELT1'] = np.abs(channelWidths[0,0]/1e6)
	header['CRVAL1'] = channelFreqs[0,0]/1e6
	header['CTYPE1'] = ('Frequency')
	header['CUNIT1'] = ('MHz')
	header['CRPIX2'] = 1
	header['CRVAL2'] = 1
	header['CDELT2'] = 1
	header['CTYPE2'] =  ('Baseline')
	header['POLAR'] =  ('XX')
	header['BUNIT'] = ('% > '+str(5)+'*rms')
	header['BTYPE'] = ('intensity')

	#write file
	hdulist.writeto('/home/maccagni/Projects/ApertiF/RFI/data/prova.fits',overwrite=True)
	hdulist.close()

	return 0

def plot_rfi_im(fitsfile):
	'''
	Plots the .fits image output of rfi_im jpg format.
	'''        
	#check if image exists
	#plot image
	fig = aplpy.FITSFigure(fitsfile,figsize=(12,8))

	#plot colorscale & colorbar
	fig.show_colorscale(aspect='auto', cmap='nipy_spectral_r',vmin=0,vmax=100)
	fig.show_colorbar()
	fig.colorbar.set_width(0.2)
	fig.colorbar.set_font(size=20, weight='medium', \
						  stretch='normal', family='sans-serif', \
						  style='normal', variant='normal')
	fig.colorbar.set_axis_label_font(size=20)
	fig.colorbar.set_axis_label_text(r'$\% > 5 \times$ r.m.s.')

	#set axis
	fig.axis_labels.set_font(size=20, weight='medium', \
							 stretch='normal', family='sans-serif', \
							 style='normal', variant='normal')
	fig.tick_labels.set_font(size=20, weight='medium', \
							 stretch='normal', family='sans-serif', \
							 style='normal', variant='normal') 
	#titleplot = self.target+': '+self.aperfi_startime+' - '+self.aperfi_endtime
	fig.savefig('/home/maccagni/Projects/ApertiF/RFI/data/prova.jpg',format='jpg')

	return 0


print ''
print '--- Working on file {0:s} ---'.format(MS[0])

t=tables.table(MS[0])
fieldIDs=t.getcol('FIELD_ID')
antennas = tables.table(MS[0]+'/ANTENNA')
ant_pos = np.array(antennas.getcol('POSITION'))
ant_wsrtnames = np.array(antennas.getcol('NAME'))
ant_names = np.arange(0,ant_wsrtnames.shape[0],1)
nant = len(ant_names)
print ant_names
#sort baselines by length
baselines = []
for i in range(0,nant-1) :
	for j in range(i+1,nant) :

		# Exclude anticorrelations
		if ant_names[i]!=ant_names[j] and ant_names[i] != 5 and ant_names[j] != 5 and ant_names[i] != 2 and ant_names[j] != 2 and ant_names[i] != 3 and ant_names[j] != 3 and ant_names[i] != 4 and ant_names[j] != 4 and ant_names[i] != 6 and ant_names[j] != 6 :
			#distances are in units of meters
			xdiff = (ant_pos[i,0]-ant_pos[j,0])*(ant_pos[i,0]-ant_pos[j,0])
			ydiff = (ant_pos[i,1]-ant_pos[j,1])*(ant_pos[i,1]-ant_pos[j,1])
			zdiff = (ant_pos[i,2]-ant_pos[j,2])*(ant_pos[i,2]-ant_pos[j,2])

			dist = np.sqrt(xdiff+ydiff+zdiff)
			baselines.append(([ant_names[i],ant_names[j]],dist))

baselines_sort = sorted(baselines, key=lambda baselines: baselines[1])  
print baselines_sort
# Define matrix of indecese of baselines                                         
blMatrix=np.zeros((nant,nant),dtype=int)
for i in range(0,len(baselines_sort)) :

	ant1 = baselines_sort[i][0][0]
	ant2 = baselines_sort[i][0][1]
	blMatrix[ant1,ant2] = i
	blMatrix[ant2,ant1] = i
print blMatrix

#load visibilities

vis = t.getcol('DATA')
print vis.shape
flag = t.getcol('FLAG')
ant1=t.getcol('ANTENNA1')
ant2=t.getcol('ANTENNA2')
#datacube visibilities sorted by baseline length, frequency, time
datacube = np.zeros([len(baselines_sort),vis.shape[1],vis.shape[0]/(len(baselines_sort))])
baseline_counter = np.zeros((nant,nant),dtype=int)
#flag unused polarizations
flag[:,:,1] = True
flag[:,:,2] = True
flag[:,:,3] = True

for i in xrange(0,vis.shape[0]):

	if (ant1[i]==0 and ant2[i] == 1) or (ant1[i] == 1 and ant2[i] == 0):
		indice=blMatrix[ant1[i],ant2[i]]
		counter=baseline_counter[ant1[i],ant2[i]]
		datacube[indice,:,counter]=np.abs(vis[i,:,0])
		baseline_counter[ant1[i],ant2[i]]+=1
	else:
		flag[i,:,0] = True

spw=tables.table(MS[0]+'/SPECTRAL_WINDOW')
channelWidths=spw.getcol('CHAN_WIDTH')
print channelWidths.shape
channelFreqs=spw.getcol('CHAN_FREQ')
print channelFreqs.shape

freq_start = channelFreqs[0]
freq_end = channelFreqs[-1]

rfifree_min = 1.180e9
rfifree_max = 1.200e9
rmsclip = 2.
rms = np.zeros([datacube.shape[0],datacube.shape[1]])
mean_array = np.zeros([datacube.shape[0],datacube.shape[1]])
flag_lim_array= np.zeros([datacube.shape[0]])

chan_min = np.argmin(np.abs(channelFreqs - rfifree_min))
chan_max = np.argmin(np.abs(channelFreqs - rfifree_max))
time_ax_len = int(datacube.shape[2])

for i in xrange(0,datacube.shape[0]):

	tmp_rms = np.nanmedian(datacube[i, chan_min:chan_max, :])
	med2 = abs(datacube[i, chan_min:chan_max, :] - tmp_rms)
	madfm = np.ma.median(med2) / 0.6744888
	flag_lim = rmsclip*madfm  
	flag_lim_array[i] = flag_lim      
	for j in xrange(0,datacube.shape[1]):

		tmpar = datacube[i,j,:]
		mean  = np.nanmean(tmpar)
		tmpar = tmpar-mean
		tmpar = abs(tmpar)
		mean_array[i,j] = mean
		#change masked values to very high number
		inds = np.where(np.isnan(tmpar))
		tmpar[inds]=np.inf
		tmpar.sort()
		index_rms = np.argmin(np.abs(tmpar - flag_lim))
		tmp_over = len(tmpar[index_rms:-1])+1
		if tmp_over == 1. :
			tmp_over = 0.
		rms[i,j] = 100.*tmp_over/time_ax_len

fitsfile = '/home/maccagni/Projects/ApertiF/RFI/data/prova.fits'

#write_freq_base(rms,channelWidths,channelFreqs,fitsfile)
#plot_rfi_im(fitsfile)


# # #Mask RFI
# for i in xrange(0,vis.shape[0]):
	
# 	if ant1[i] == ant2[i] or ant1[i] == 5 or ant2[i] == 5:
# 		continue
# 	else:
# 		indice=blMatrix[ant1[i],ant2[i]]
# 		flag_clip = flag_lim_array[indice]
	
# 	for j in xrange(0,vis.shape[1]):
# 		mean = mean_array[indice,i]
# 		if (np.abs(vis[i,j,0]) - mean ) > flag_clip:
			
# 			flag[i,j,0] = True


newtable = 	'/home/maccagni/prova_2ant.MS'
tout = t.copy('/home/maccagni/prova_2ant.MS')       # make deep copy
tout.close()
t.close()
tout = tables.table(newtable, readonly=False)
#coldmi = tout.getdminfo('DATA')     # get dminfo of existing column
#coldmi["NAME"] = 'CORRECTED_DATA'               # give it a unique name
#tout.addcols(tables.maketabdesc(tables.makearrcoldesc("CORRECTED_DATA",0., ndim=2)),
#           coldmi)
#tout.putcol('CORRECTED_DATA',vis)
tout.putcol('FLAG',flag)
tout.close()


print 'NORMAL TERMINATION'