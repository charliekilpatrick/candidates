import warnings
warnings.filterwarnings('ignore')
from astropy.io import ascii
from astropy.table import Table, Column, vstack, unique, hstack
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck15 as cosmo
from astropy import units as u
from bs4 import BeautifulSoup
from scipy.optimize import curve_fit
import copy, sys, requests, os, time, pandas, progressbar
import healpy as hp
import numpy as np
import astropy_healpix as ah
from ligo.skymap import distance as lvc_distance

def estimate_area(probability, event):
            gw_map = f'data/{event}/{event}_skymap.fits.gz'
            map_test = Table.read(gw_map)

            level, ipix = ah.uniq_to_level_ipix(map_test['UNIQ'])
            nside = ah.level_to_nside(level)
            area_per_pix = hp.pixelfunc.nside2pixarea(nside, degrees=True)

            hpx = map_test['PROBDENSITY']

            sorted_pixels = np.flip(np.argsort(hpx))

            dist_map = map_test['DISTMU']
            sigma_map = map_test['DISTSIGMA']
            header = map_test.meta

            total_prob = 0.0
            total_area = 0.0
            for i,pix in enumerate(sorted_pixels):
                if total_prob > probability:
                    break
                else:
                    total_prob += hpx[pix] * area_per_pix[pix] * (np.pi/180.0)**2
                    total_area += area_per_pix[pix]

            m = 'Total area for {prob} probability is {area} deg^2'
            print(m.format(prob=probability, area=total_area))

estimate_area(0.90)
estimate_area(0.95)
estimate_area(0.99)
estimate_area(0.995)
