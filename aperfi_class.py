#!/usr/bin/env python
__author__ = "Filippo Maccagni, Tom Osterloo"
__copyright__ = "ASTRON"
__email__ = "maccagni@astron.nl"

# Import modules
import ConfigParser
import logging
from libs import lib
import os
import aplpy
import  astropy.io.fits as pyfits
import string
from matplotlib import gridspec
from matplotlib import pyplot as plt
#import convert
#import scal
import setinit_an
import numpy as np
import pyrap.tables as tables
import numpy.ma as ma
import aplpy
from astropy.io import fits
#from apercal import convert


class aperfi_class:
    '''
    
    Class to investigate the RFI behaviour during observations

    '''

    def __init__(self, file=None, **kwargs):
        '''
    
        Set logger for spectrum extraction
        Find config file
        If not specified by user load default.cfga
    
        '''
        
        self.logger = logging.getLogger('RFI')
        config = ConfigParser.ConfigParser() # Initialise the config parser
        if file != None:
            config.readfp(open(file))
            self.logger.info('### Configuration file ' + file + ' successfully read! ###')
        else:
            config.readfp(open(os.path.realpath(__file__).rstrip('aperfi_class.pyc') + 'default.cfga'))
            self.logger.info('### No configuration file given or file not found! Using default values! ###')
        for s in config.sections():
            for o in config.items(s):
                setattr(self, o[0], eval(o[1]))
        self.default = config # Save the loaded config file as defaults for later usage
        setinit_an.setinitdirs_rfi(self)
        setinit_an.setinitfiles_rfi(self)

    def go(self):
        '''

        Executes the whole spectrum extraction process as follows:
        1: load_from_ms
        2: baselines_from_ms 
        3: priors_flag
        4: rfi_flag

        '''

        self.logger.info("########## STARTING RFI analysis ##########")
        self.load_from_ms()
        self.baselines_from_ms()
        self.logger.info("---- Baselines written -----")
        self.priors_flag()
        self.logger.info("---- Bad antennas and autocorrelations flagged ----")
        self.rfi_flag()
        self.logger.info("---- RFI found ----")
        self.logger.info("########## END RFI ANALYSIS ##########")

    def load_from_ms(self):
        '''

        Loads important columns from MS file
        From MS: 
        Field_ID, Data, Flag, Antenna1, Antenna2,
        From MS/ANTENNA:
        Position,Name
        From MS/SPECTRAL_WINDOW
        Chan_width, Chan_freq

        '''


        t=tables.table(self.msfile)

        self.fieldIDs=t.getcol('FIELD_ID')
        self.vis = t.getcol('DATA')
        print self.vis.shape
        self.flag = t.getcol('FLAG')
        self.ant1=t.getcol('ANTENNA1')
        self.ant2=t.getcol('ANTENNA2')
        
        t.close()

        #
        antennas = tables.table(self.msfile +'/ANTENNA')
        self.ant_pos = np.array(antennas.getcol('POSITION'))
        self.ant_wsrtnames = np.array(antennas.getcol('NAME'))
        self.ant_names = np.arange(0,self.ant_wsrtnames.shape[0],1)
        self.nant = len(self.ant_names)
        self.logger.info("\tTotal number of antennas\t:"+str(self.nant))
        self.logger.info("\tAntenna names\t"+str(self.ant_names))

        print "\tTotal number of antennas\t:"+str(self.nant)
        print "\tAntenna names\t"+str(self.ant_names)
        

        antennas.close()

        spw=tables.table(self.msfile+'/SPECTRAL_WINDOW')
        self.channelWidths=spw.getcol('CHAN_WIDTH')
        self.channelFreqs=spw.getcol('CHAN_FREQ')
        self.logger.info("\tBandwidth\t:"+str(self.channelWidths))
        self.logger.info("\tStart Frequency\t:"+str(self.channelFreqs[0][0]))
        self.logger.info("\tEnd Frequency\t:"+str(self.channelFreqs[-1][0]))

        print "\tBandwidth [kHz]\t:"+str(self.channelWidths[0][0]/1e3)
        print "\tStart Frequency [GHz]\t:"+str(self.channelFreqs[0][0]/1e9)
        print "\tEnd Frequency [GHz]\t:"+str(self.channelFreqs[-1][-1]/1e9)
        spw.close()



    def baselines_from_ms(self):
        '''

        Reads which baselines were used in the observations
        Stores them sorted by lenght in the array baselines_sort
        Creates the Matrix to analize visibilites on each baseline separately

        '''

        #sort baselines by length
        baselines = []
        for i in range(0,self.nant-1) :
            for j in range(i+1,self.nant) :

                # Exclude anticorrelations
                if self.ant_names[i]!=self.ant_names[j] and any(x != self.ant_names[i] for x in self.aperfi_badant) and any(x != self.ant_names[i] for x in self.aperfi_badant):
                    #distances are in units of meters
                    xdiff = (self.ant_pos[i,0]-self.ant_pos[j,0])*(self.ant_pos[i,0]-self.ant_pos[j,0])
                    ydiff = (self.ant_pos[i,1]-self.ant_pos[j,1])*(self.ant_pos[i,1]-self.ant_pos[j,1])
                    zdiff = (self.ant_pos[i,2]-self.ant_pos[j,2])*(self.ant_pos[i,2]-self.ant_pos[j,2])

                    dist = np.sqrt(xdiff+ydiff+zdiff)
                    baselines.append(([self.ant_names[i],self.ant_names[j]],dist))

        self.baselines_sort = sorted(baselines, key=lambda baselines: baselines[1])  

        # Define matrix of indecese of baselines                                         
        self.blMatrix=np.zeros((self.nant,self.nant),dtype=int)
        for i in range(0,len(self.baselines_sort)) :

            ant1 = self.baselines_sort[i][0][0]
            ant2 = self.baselines_sort[i][0][1]
            self.blMatrix[ant1,ant2] = i
            self.blMatrix[ant2,ant1] = i

        print self.blMatrix

    def priors_flag(self):
        '''

        Flags YY,XY,YX polarizations
        Flags autocorrelations
        Flags bad antennas set by aperfi_badant = [ x, y, z ]
        Stores them sorted by lenght in the array baselines_sort

        '''

        self.datacube = np.zeros([len(self.baselines_sort),self.vis.shape[1],self.vis.shape[0]/(len(self.baselines_sort))])
        baseline_counter = np.zeros((self.nant,self.nant),dtype=int)

        #flag unused polarizations
        self.flag[:,:,1] = True #YY
        self.flag[:,:,2] = True #XY
        self.flag[:,:,3] = True #YX

        #flag autocorrelations and bad antennas
        for i in xrange(0,self.vis.shape[0]):

            if self.ant1[i] == self.ant2[i]:
                self.flag[i,:,0] = True
                
            elif (any(x == self.ant1[i] for x in self.aperfi_badant) or any(x == self.ant2[i] for x in self.aperfi_badant)):
                self.flag[i,:,0] = True

            else:
                a1 = self.ant1[i]
                a2 = self.ant2[i]
                indice=self.blMatrix[a1,a2]
                # Count how many visibilities per baseline
                counter=baseline_counter[a1,a2]
                # Put amplitude of visibility
                # In the right place in the new array
                self.datacube[indice,:,counter]=np.abs(self.vis[i,:,0])
                # Update the number of visibility in that baseline
                baseline_counter[a1,a2]+=1

        #self.datacube = np.transpose(self.datacube, (1, 0, 2)) 


    def find_rfi(self):
        '''

        For each baseline finds all signal above rms*aperfi_clip
        Creates a cube of visibilities sorted by baseline_lengh , frequency, time
        Stores them sorted by lenght in the array baselines_sort
        Creates the Matrix to analize visibilites on each baseline separately

        '''


        self.rms = np.zeros([self.datacube.shape[0],self.datacube.shape[1]])
        self.mean_array = np.zeros([self.datacube.shape[0],self.datacube.shape[1]])
        self.flag_lim_array= np.zeros([self.datacube.shape[0]])

        chan_min = np.argmin(np.abs(self.channelFreqs[0] - self.aperfi_rfifree_min))
        chan_max = np.argmin(np.abs(self.channelFreqs[0] - self.aperfi_rfifree_max))
        time_ax_len = int(self.datacube.shape[2])
        for i in xrange(0,self.datacube.shape[0]):
            tmp_rms = np.nanmedian(self.datacube[i, chan_min:chan_max, 0])
            med2 = abs(self.datacube[i, chan_min:chan_max, 0] - tmp_rms)
            madfm = np.ma.median(med2) / 0.6744888
            flag_lim = self.aperfi_rmsclip*madfm  
            self.flag_lim_array[i] = flag_lim    

            for j in xrange(0,self.datacube.shape[1]):
                tmpar = self.datacube[i,j,:]
                mean  = np.nanmean(tmpar)
                tmpar = tmpar-mean
                tmpar = abs(tmpar)
                self.mean_array[i,j] = mean
                #change masked values to very high number
                #inds = np.where(np.isnan(tmpar))
                tmpar[np.isnan(tmpar)]=np.inf
                tmpar.sort()
                index_rms = np.argmin(np.abs(tmpar - flag_lim))
                tmp_over = len(tmpar[index_rms:-1])+1
                if tmp_over == 1. :
                    tmp_over = 0.
                self.rms[i,j] = 100.*tmp_over/time_ax_len

        self.write_freq_base()

        return 0


    def rfi_flag(self):
        '''
        Creates a new MS file where RFI has been flagged in the FLAG column
        '''

        #self.flag_lim_array = np.zeros([self.vis.shape[0]])
        for i in xrange(0,self.vis.shape[0]):
            
            if any(x == self.ant1[i] for x in self.aperfi_badant) or any(x == self.ant2[i] for x in self.aperfi_badant):
                continue
            else:
                indice=self.blMatrix[self.ant1[i],self.ant2[i]]
                flag_clip = self.flag_lim_array[indice]
            
            for j in xrange(0,self.vis.shape[1]):
                mean = self.mean_array[indice,i]
                if (np.abs(self.vis[i,j,0]) - mean ) > flag_clip:
                
                    self.flag[i,j,0] = True
        
        t=tables.table(self.msfile)
        tout = t.copy(self.rfifile)       # make deep copy
        tout.close()
        t.close()
        tout = tables.table(self.rfifile, readonly=False)
        #coldmi = tout.getdminfo('DATA')     # get dminfo of existing column
        #coldmi["NAME"] = 'CORRECTED_DATA'               # give it a unique name
        #tout.addcols(tables.maketabdesc(tables.makearrcoldesc("CORRECTED_DATA",0., ndim=2)),
        #           coldmi)
        #tout.putcol('CORRECTED_DATA',vis)
        tout.putcol('FLAG',self.flag)
        tout.close()

    def make_psf(self) :
        '''
        use wsclean to predict the psf of the observation
        '''

        command1 = '''wsclean -name psfonly -mem 100 -no-dirty -weight natural -super-weight 1.0'''
        command2 = '''-weighting-rank-filter-size 16 -size 512 512 -scale 3.0asec -channels-out 1'''
        command3 = ''' -grid-mode kb -kernel-size 7 -oversampling 63 -make-psf-only -pol xx -intervals-out 1'''
        command4 = ''' -data-column DATA -gain 0.1 -mgain 1.0 -multiscale-scale-bias 0.6 -fit-beam -no-fit-beam '''+self.rfifile

        command = command1+command2+command3+command4

        os.system(command)


    def write_freq_base(self) :
        '''
        
        Writes an image.fits of visibilities ordered by frequency, baseline.
        Baselines are ordered by their length.
        
        '''
        #reverse frequencies if going from high-to-low         
        
        #set fits file
        hdu = fits.PrimaryHDU(self.rms)
        hdulist = fits.HDUList([hdu])
        header = hdulist[0].header
        #write header keywords               
        header['CRPIX1'] = 1
        header['CDELT1'] = np.abs(self.channelWidths[0,0]/1e6)
        header['CRVAL1'] = self.channelFreqs[0,0]/1e6
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
        fits.writeto(self.rfi_freq_base,self.rms,header)
        hdulist.close()

        return 0

    def rfi_frequency(self):
        '''
        Determines the rfi per frequency channel. Saves results in table rfi_table.fits
        For each channel the flag % and factor of noise increase are stored for all, long and short baselines
        Long and short baselines are separated in half, depending on the number of baselines
        '''

        #open file
        if os.path.exists(self.rfi_freq_base) == False:
            self.logger.error('### Image of RFI sorted by frequency over baseline lenght does not exist ###')    
            self.logger.error('### Run aperfi.rfi_im() first ###')  
        else:    
            
            # read data and header
            hdulist = pyfits.open(self.rfi_freq_base)  # read input                
            datacube = hdulist[0].data    
            prihdr = hdulist[0].header

            #set array of frequencies
            self.freqs = (np.linspace(1, datacube.shape[1], datacube.shape[1])\
                         - prihdr['CRPIX1'])*prihdr['CDELT1'] + prihdr['CRVAL1']
            
            # set y-array
            rms_lin = np.zeros([datacube.shape[1]])    
            flag_lin = np.zeros([datacube.shape[1]])    
            rms_lin_long = np.zeros([datacube.shape[1]]) + np.sqrt(2.)          
            rms_lin_short = np.zeros([datacube.shape[1]]) + np.sqrt(2.)   
            flag_lin_long = np.zeros([datacube.shape[1]]) + 50.          
            flag_lin_short = np.zeros([datacube.shape[1]]) + 50.

            for i in xrange(0,datacube.shape[1]):
                
                flag_lin_tmp = np.divide(np.sum(datacube[:,i]),datacube.shape[0])
                flag_lin[i] = flag_lin_tmp

                shortbase=datacube[:int(datacube.shape[0]/2),i]
                longbase = datacube[int(datacube.shape[0]/2):,i]               
                
                rms_lin_tmp = 1.-np.divide(np.divide(np.sum(datacube[:,i]),datacube.shape[0]),100.)
                rms_lin[i] = np.divide(1.,np.sqrt(rms_lin_tmp))

                flag_lin_tmp = np.divide(np.sum(shortbase),len(shortbase))
                flag_lin_short[i] = flag_lin_tmp
                rms_lin_tmp_short = 1.-np.divide(np.divide(np.sum(shortbase),len(shortbase)),100.)
                rms_lin_short[i] *= np.divide(1.,np.sqrt(rms_lin_tmp_short))

                flag_lin_tmp = np.divide(np.sum(longbase),len(longbase))
                flag_lin_long[i] = flag_lin_tmp
                rms_lin_tmp_long = 1.-np.divide(np.divide(np.sum(longbase),len(longbase)),100.)
                rms_lin_long[i] *= np.divide(1.,np.sqrt(rms_lin_tmp_long))
        
            # save fits table        
            c1 = pyfits.Column(name='frequency', format='D', unit='MHz', array=self.freqs)
            c2 = pyfits.Column(name='flag', format='D', unit='-', array=flag_lin)
            c3 = pyfits.Column(name='noise_factor', format='D', unit = '-', array=rms_lin)
            c4 = pyfits.Column(name='flag_short', format='D', unit='-', array=flag_lin_short)
            c5 = pyfits.Column(name='noise_factor_short', format='D', unit = '-', array=rms_lin_short)
            c6 = pyfits.Column(name='flag_long', format='D', unit='-', array=flag_lin_long)
            c7 = pyfits.Column(name='noise_factor_long', format='D', array=rms_lin_long)        

            fits_table = pyfits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7])    
            
            fits_table.writeto(self.rfi_table, overwrite = True)
            

    def plot_rfi_im(self):      
        '''
        
        Plots the .fits image output of rfi_im jpg format.
        
        '''        
        
        #check if image exists
        #plot image
        fig = aplpy.FITSFigure(self.rfi_freq_base,figsize=(12,8))

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
        plt.savefig(self.rfi_freq_base_plot,format='png' ,overwrite=True)


    def plot_noise_frequency(self):
        '''
        Plots the noise or % of rfi per frequency channel for all, long and short baselines.
        In default.cfga
            aperfi_noise = 'rfi' (or flag)
            aperfi_plot_long_short = False
        '''

        #open file
        if os.path.exists(self.rfi_table) == False:
            self.logger.error('### Table of RFI and flags of visibilities does not exist ###')    
            self.logger.error('### Run aperfi.rfi_frequency() first ###')  
        else:  


            t = pyfits.open(self.rfi_table)
            data_vec = t[1].data
            cols = t[1].columns
            
            freqs = np.array(data_vec['frequency'],dtype=float)
            flags = np.array(data_vec['flag'],dtype=float)
            noise_factor = np.array(data_vec['noise_factor'],dtype=float)
            noise_factor_long = np.array(data_vec['noise_factor_long'],dtype=float)
            flags_long = np.array(data_vec['flag_long'],dtype=float)
            noise_factor_short = np.array(data_vec['noise_factor_short'],dtype=float)
            flags_short = np.array(data_vec['flag_short'],dtype=float)

           
            #if self.aperfi_noise == 'noise':
            #    self.predicted_noise_channel()
            #    noise_all = noise_factor*self.noise_freq
            #    noise_short = noise_factor_short*self.noise_freq
            #    noise_long = noise_factor_long*self.noise_freq
            if self.aperfi_plot_noise == 'rfi':
                noise_all = noise_factor
                noise_short = noise_factor_short
                noise_long = noise_factor_long          
            if self.aperfi_plot_noise == 'flag':
                noise_all = flags
                noise_long = flags_long
                noise_short = flags_short


            # initialize plotting parameters
            params = {'font.family'         :' serif',
                      'font.style'          : 'normal',
                      'font.weight'         : 'medium',
                      'font.size'           : 20.0,
                      'text.usetex': True,
                      'text.latex.unicode': True
                       }
            plt.rcParams.update(params)
            
            # initialize figure
            fig = plt.figure(figsize =(14,8))
            fig.subplots_adjust(hspace=0.0)
            gs = gridspec.GridSpec(1, 1)
            plt.rc('xtick', labelsize=20)

            # Initialize subplots
            ax1 = fig.add_subplot(gs[0])
            ax1.set_xlabel(r'Frequency [MHz]',fontsize=20)
            
            if self.aperfi_noise != 'flag':
                ax1.set_yscale('log', basey=10)
        
            #define title output                     
            
            #plot
            label_all = 'All baselines' 
            label_long = 'Long baselines' 
            label_short = 'Short baselines' 

            if self.aperfi_plot_long_short == True:
                ax1.step(freqs,noise_short, where= 'pre', color='red', linestyle='-',label=label_short)
                ax1.step(freqs,noise_long, where= 'pre', color='blue', linestyle='-',label=label_long)
                out_plot = out_plot+'_sl_'

            ax1.step(freqs,noise_all, where= 'pre', color='black', linestyle='-',label=label_all)

            #titleplot = self.target+': '+self.aperfi_startime+' - '+self.aperfi_endtime
            #ax1.set_title(titleplot)
            
            # set axis, legend ticks

            ax1.set_xlim([np.min(freqs)-5,np.max(freqs)+5])
            xticks_num = np.linspace(int(self.channelFreqs[0]),int(self.channelFreqs[-1]),10,dtype=int)
            ax1.set_xticks(xticks_num)

            if self.aperfi_plot_noise == 'rfi':
                ax1.set_yticks([1,round(np.sqrt(2),2),2,3,5,10,50]) 
                ax1.set_ylabel(r'Factor of noise increase')
                ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

            # if self.aperfi_noise == 'noise':
            #     ax1.set_yticks([1,2,3,5,10,50]) 
            #     ax1.set_ylabel(r'Predicted noise [mJy beam$^{-1}$]')     
            #     out_plot = out_plot+'_noise'+self.aperfi_plot_format    
            #     ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

            if self.aperfi_plot_noise == 'flag':
                ax1.set_ylabel(r'$\% >$ '+str(self.aperfi_rmsclip)+'*rms') 
            
            legend = plt.legend()
            legend.get_frame().set_edgecolor('black')

            # Save figure to file
            plt.savefig(self.rfi_freq_plot,format='png',overwrite = True)
                            

    #######################################################################
    ##### Utilities                                                   #####
    #######################################################################
    
    def show(self, showall=False):
        '''
        
        show: Prints the current settings of the pipeline. 
              Only shows keywords, which are in the default analysis config file apercal/third_party/default.cfga
        
        showall=True : see all current settings of default.cfga instead of only the ones from the current class
        
        '''
        setinit_an.setinitdirs_abs(self)
        config = ConfigParser.ConfigParser()
        config.readfp(open(self.apercaldir + '/third_party/default.cfga'))
        for s in config.sections():
            if showall:
                print(s)
                o = config.options(s)
                for o in config.items(s):
                    try:
                        print('\t' + str(o[0]) + ' = ' + str(self.__dict__.__getitem__(o[0])))
                    except KeyError:
                        pass
            else:
                if s == 'RFI':
                    print(s)
                    o = config.options(s)
                    for o in config.items(s):
                        try:
                            print('\t' + str(o[0]) + ' = ' + str(self.__dict__.__getitem__(o[0])))
                        except KeyError:
                            pass
                else:
                    pass
    
    
    def reset(self):
        '''
        
        Resets the current step and remove all generated data. 
        Be careful! Deletes all data generated in this step!
        
        '''
        self.logger.warning('### Deleting all preflagged data. You might need to copy over the raw data to the raw subdirectory again. ###')
        self.director('ch', self.rawdir)
        self.director('rm', self.rawdir + '/*')



