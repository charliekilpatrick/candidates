import os
import sys
import shutil
from astropy.utils.data import download_file

files=['hlsp_ps1-strm_ps1_gpc1_m33-m28_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_m28-m25_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_m25-m22_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_m22-m18_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_m18-m15_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_m15-m12_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_m12-m09_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_m09-m05_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_m05-m02_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_m02-p01_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p01-p05_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p05-p08_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p08-p11_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p11-p15_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p15-p18_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p18-p21_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p21-p25_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p25-p28_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p28-p31_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p31-p34_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p34-p38_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p38-p41_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p41-p44_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p44-p48_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p48-p51_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p51-p54_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p54-p58_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p58-p61_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p61-p64_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p64-p69_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p69-p77_multi_v1_cat.csv.gz',
'hlsp_ps1-strm_ps1_gpc1_p77-p90_multi_v1_cat.csv.gz']

dl_dir = ''
if len(sys.argv)>1:
    dl_dir = sys.argv[1]
else:
    dl_dir = os.getcwd()

base_url = 'http://archive.stsci.edu/hlsps/ps1-strm/'

for file in files:

    url = base_url + file
    outfile = os.path.join(dl_dir, file)

    if not os.path.exists(outfile):
        print(f'Getting: {outfile}')
        dat = download_file(url, cache=False, show_progress=True, timeout=120)
        shutil.move(dat, outfile)
