#! /home/apercal/pipeline/bin/python python
# Import MoD_AbS modules
import sys, os

sys.path.append('/home/maccagni/programs/Apertif_modules/')
from RFInder import RFInder as RFI



filename = sys.argv[-1]
rfinder = RFI(filename)
rfinder.show()


rfinder.load_from_ms()
rfinder.baselines_from_ms()
rfinder.priors_flag()
rfinder.find_rfi()
rfinder.rfi_flag()
rfinder.plot_rfi_im()

#convert flagged dataset

#line = "exportuvfits(vis='"+aperfi.rfifile+"', fitsfile='"+aperfi.rfifile.rstrip('MS')+"UVFITS')"
#f=open('export_script.py','w')
#f.write(line)

#print rfinder.rfifile
#command = '''wsclean -name psfonly -mem 100 -no-dirty -weight natural -super-weight 1.0 -weighting-rank-filter-size 16 -size 1024 1024 -scale 2.6asec -channels-out 1 -grid-mode kb -kernel-size 7 -oversampling 63 -make-psf-only -pol xx -intervals-out 1 -data-column DATA -gain 0.1 -mgain 1.0 -multiscale-scale-bias 0.6 -fit-beam -no-fit-beam '''+aperfi.rfifile

#os.system(command)
#command = ''' casa -c 'export_script.py' '''
#os.system(command)


#ms = casac.casac.ms()
#ms.open(aperfi.rfifile)
#ms.tofits(aperfi.rfifile.rstrip('MS') + 'UVFITS', column='DATA', combinespw=True, padwithflags=True, multisource=True, writestation=True)

