import requests
import numpy as np
import os
import sys
import shutil
import re
from astropy.utils.data import download_file

base_urls=['https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9/south/sweep/9.0-photo-z/',
           'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/',
           'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9/north/sweep/9.0-photo-z/'
           'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9/north/sweep/9.0/']

dl_dir = ''
if len(sys.argv)>1:
    dl_dir = sys.argv[1]
else:
    dl_dir = os.getcwd()

for url in base_urls:

    r = requests.get(url)

    if r.status_code==200:

        data = re.findall('sweep\-.*?\.fits',r.text)
        data = np.unique(data)

        for d in data:
                outfile = os.path.join(dl_dir, d)
                dl_url = os.path.join(url, d)

                if not os.path.exists(outfile):

                    print(f'Getting: {outfile}')
                    dat = download_file(dl_url, cache=False,
                        show_progress=True, timeout=120)
                    shutil.move(dat, outfile)



