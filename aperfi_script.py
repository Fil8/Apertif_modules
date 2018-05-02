#! /home/apercal/pipeline/bin/python python

from third_party.aperfi_class import aperfi_class as aperfi
import sys, os

filename = sys.argv[-1]
aperfi = aperfi(filename)
aperfi.show()

print aperfi.rfidir

aperfi.load_from_ms()
aperfi.baselines_from_ms()
aperfi.priors_flag()
aperfi.rfi_flag()

#convert flagged dataset

#line = "exportuvfits(vis='"+aperfi.rfifile+"', fitsfile='"+aperfi.rfifile.rstrip('MS')+"UVFITS')"
#f=open('export_script.py','w')
#f.write(line)

print aperfi.rfifile
command = '''wsclean -name psfonly -mem 100 -no-dirty -weight natural -super-weight 1.0 -weighting-rank-filter-size 16 -size 1024 1024 -scale 2.6asec -channels-out 1 -grid-mode kb -kernel-size 7 -oversampling 63 -make-psf-only -pol xx -intervals-out 1 -data-column DATA -gain 0.1 -mgain 1.0 -multiscale-scale-bias 0.6 -fit-beam -no-fit-beam '''+aperfi.rfifile

os.system(command)
#command = ''' casa -c 'export_script.py' '''
#os.system(command)


#ms = casac.casac.ms()
#ms.open(aperfi.rfifile)
#ms.tofits(aperfi.rfifile.rstrip('MS') + 'UVFITS', column='DATA', combinespw=True, padwithflags=True, multisource=True, writestation=True)

