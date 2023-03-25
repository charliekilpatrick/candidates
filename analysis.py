import warnings, sys, os
warnings.filterwarnings('ignore')
from astropy.io import ascii
from astropy.table import Table, Column, vstack
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u

from astropy.utils.data import download_file
import healpy as hp # for working with HEALPix files

from utilities.utilities import *

event = sys.argv[1]
yse_table = sys.argv[2]
redo=False

shibboleth = '/home/ckilpatrick/scripts/shibboleth'

# Constraints on candidates and search parameters
constraints = {
    'time': [0, 14] * u.day,     # Range in days for viable candidates
    'redshift_range': [0.0001, 0.35],
    'probability': 0.95,         # Percentile range for viable candidates
    'distance_percentile': 0.95,
    'search_radius': {
        'mpc': 20.0 * u.arcsec,
        'gaia': 1.0 * u.arcsec,
        'asassn': 1.0 * u.arcsec,
        'ps1dr2': 1.0 * u.arcsec,
        'galaxy': 30.0 * u.arcsec,
        'galaxy_proj': 3.0e5 * u.parsec,
    },
    'photometry': {
        'bright': -18.0,
        'decline': 0.1,
        'color': -1.0,
    }
}

# Table names
tbl = {
    'candidate': 'initial_candidate_table',
    'internal': 'internal_candidate_table',
    'classification': 'classification_table',
    'spec': 'spec_z_table',
    'ps1dr2': 'ps1dr2_table',
    'asassn': 'asassn_table',
    'ps1strm': 'ps1strm_table',
    '2mpz': '2mpz_table',
    'legacy': 'legacy_table',
}

# Directory structure and organization
event_dir = os.path.join('data', event) ; tmp = os.path.join('data', event, 'tmp')
if not os.path.isdir(event_dir):
    try:
        os.makedirs(event_dir)
    except:
        error = 'ERROR: could not create directory {dir}.  Exiting...'
        print(error.format(event_dir))
        sys.exit()
for key in tbl.keys():
    tbl[key] = os.path.join(event_dir, tbl[key])

url=f'https://gracedb.ligo.org/api/superevents/{event}/files/LALInference.fits.gz'
event_file = os.path.join(event_dir, event+'_LALInference.pkl')
if os.path.exists(event_file):
    data = pickle.load(event_file)
    header = data['header']
    prob = data['prob']
    table = data['table']
else:
    filename = download_file(url, cache=True)
    prob, header = hp.read_map(filename, h=True)
    header = dict(header)
    table = Table.read(url)

    data = {'header': header, 'prob': prob, 'table': table}
    pickle.dumps(event_file)

# Storing event and analysis specific parameters
kwargs = {'alert_time': Time(header['DATE-OBS']),
          'meta': [],
          'reference': None,
          'reference_name': '',
          'yse': yse_table,
          'gw_map_url': url,
          'gw_map': prob,
          'gw_map_table': table,
          'gw_map_header': header,
          'redo': redo,
          'avoid_filters': ['G','B','V','U','UVM2','UVW2',
                            'UVW1','orange','cyan','R'],
          'redshift': {
            'methods': ['spec','ps1strm','2mpz','legacy'],
            #'methods': ['spec','legacy'],
            'use_name': 'use_redshift',
           },
           'photyse': 106,
           # Analysis procs that don't require crossmatching
           'no_crossmatch': ['import candidates', 'time', 'prob','redshift',
                'class','photometry'],
            'outtable': {
                'ra_precision': 3,
                'dec_precision': 2,
                'sort': ['name_len','name'],
            }
}
kwargs.update(tbl)
kwargs.update(constraints)

steps = [import_candidates, check_time, check_prob, check_mpc, check_gaia,
    check_class, check_ps1dr2, check_asassn, check_redshift, check_photometry]

internal = None
prob = None

internal_candidate_file = os.path.join('data', event, 'internal_candidate_table.dat')
probability_file = os.path.join('data',event,'prob_table.dat')

if os.path.exists(internal_candidate_file):
    internal=Table.read(internal_candidate_file, format='ascii')

if os.path.exists(probability_file):
    prob=Table.read(probability_file, format='ascii')

table = None
masks = []
for i,proc in enumerate(steps):
    table = add_data(table, proc, **kwargs)
    remove_rows=[]
    if i==0:
        for j,row in enumerate(table):
            if row['name'].startswith('SSS') or row['name'].startswith('TGW'):
                if internal and row['name'] not in internal['name']:
                    remove_rows.append(j)
    table.remove_rows(remove_rows)
    remove_rows=[]
    if i==0:
        for j,row in enumerate(table):
            if prob and row['name'] not in prob['name']:
                remove_rows.append(j)
    table.remove_rows(remove_rows)
    if i>0:
        masks.append(table[format_proc(proc)+'_mask'].data)
    if len(masks)>0:
        mcopy = copy.copy(masks)
        mcopy = np.array([np.array([bool(f) for f in m]) for m in mcopy])
        mask = np.any(np.array(mcopy), axis=0)
        subtable=table[~mask]

        outdatafile = os.path.join(event_dir, format_proc(proc)+'_table.dat')
        subtable.write(outdatafile, format='ascii', overwrite=True)

    step_message(i, proc, get_kwargs(table, masks=masks))

table.meta['use_masks']=masks

output_latex_table(table, **kwargs)
