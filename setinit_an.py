import os

def setinitdirs_abs(self):
    '''
    
    Creates the directory names for the subdirectories to make scripting easier
    
    '''
    
    self.linedir = self.basedir + self.beam + '/' + self.linesubdir
    self.contdir = self.basedir + self.beam + '/' + self.contsubdir
    self.absdir = self.basedir + self.beam  + '/abs/' 
    self.specdir = self.absdir+'spec/'    
    self.plotdir = self.absdir+'plot/'    
    
def setinitfiles_abs(self):
    '''
    
    Creates the file names for spectral analisys
    
    '''
    # Name the datasets -> t
    # !!!! to change according to the final products of previous pipeline
    self.data_cube = self.cubename+'.fits'
    self.data_cube_mir = self.cubename+'.mir'
    self.cont_im = self.contname+'.fits'
    self.cont_im_mir = self.contname+'.mir'
    self.src_imsad_out = 'cont_src_imsad.txt'
    self.src_list_csv = 'cont_src.csv'

def setinitdirs_rfi(self):
    '''
    
    Creates the directory names for the subdirectories to make scripting easier
    
    '''

    self.rawdir = self.basedir + self.beam + '/' + self.rawsubdir+'/'

    self.rfidir = self.basedir + self.beam + '/' + 'rfi/' 

    if os.path.exists(self.rfidir) == False:
             os.makedirs(self.rfidir)

    self.rfiplotdir = self.rfidir+'plot/' 

    if os.path.exists(self.rfiplotdir) == False:
             os.makedirs(self.rfiplotdir)


def setinitfiles_rfi(self):
    '''
    
    Creates the file names for spectral analisys
    
    '''
    
    # Name the datasets -> t
    # !!!! to change according to the final products of previous pipeline
    self.msfile = self.rawdir+self.msname[0]
    self.rfifile = self.rfidir+'rfi_flagged_vis.MS'
    self.rfi_freq_base = self.rfidir+'freq_base.fits'
    self.rfi_freq_base_plot = self.rfiplotdir+'freq_base.png'
    self.rfimsfile = self.rfidir+'rfi_flagged.MS'
    self.rfi_table = self.rfidir+'rfi_table.fits'
    self.rfi_freq_plot = self.rfiplotdir+'freq_fri.png'



