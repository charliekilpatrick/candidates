from astropy import units as u
import os
import requests
import pickle
import pathlib
import shutil

from astropy.io import ascii
from astropy.table import Table, Column, vstack
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.utils.data import download_file

import healpy as hp # for working with HEALPix files

def astcheck_dependency_check():

    astcheck_check = shutil.which('astcheck')
    mpc2sof_check = shutil.which('mpc2sof')

    if astcheck_check is None or mpc2sof_check is None:
        m='ERROR: download and make astcheck and mpc2sof from the github below '
        m+='to enable astcheck functionality for candidate analysis. \n\n\n'
        m+='https://github.com/Bill-Gray/lunar'

        raise Exception(m)


def add_constraints(args):

    # Constraints on candidates and search parameters
    constraints = {
        'time': [0, args.max_date] * u.day, # Range (days) for viable candidates
        'redshift_range': args.redshift,
        'probability': args.max_prob, # Percentile range for viable candidates
        'distance_percentile': args.distance_distribution,
        'search_radius': {
            'mpc': args.mpc_radius * u.arcsec,
            'astcheck': args.astcheck_radius * u.arcsec,
            'gaia': args.gaia_radius * u.arcsec,
            'asassn': args.asassn_radius * u.arcsec,
            'ps1dr2': args.ps1_radius * u.arcsec,
            'galaxy': args.gal_radius * u.arcsec,
            'tns': args.tns_radius * u.arcsec,
            'yse': args.yse_radius * u.arcsec,
            'galaxy_proj': args.max_gal_radius * u.parsec,
        },
        'photometry': {
            'bright': args.phot_bright,
            'decline': args.phot_min_decline,
            'color': args.phot_color,
        }
    }

    return(constraints)

def get_candidate_table(sources):

    table = None
    if os.path.exists(sources):
        try:
            table = ascii.read(sources)
        except:
            raise Exception(f'ERROR: ascii.read could not read {sources}!')
    else:
        try:
            r = requests.get(sources)
            table = ascii.read(r.text)
        except:
            raise Exception(f'ERROR: could not download or parse data from {sources}!')

    if table is not None and isinstance(table, Table):
        keys = list(table.keys())
        for key in keys:
            if key!=key.lower():
                table.rename_column(key, key.lower())

        if 'mjd' in keys and 'discovery_date' not in keys:
            times = [Time(r['mjd'], format='mjd') for r in table]
            times = [t.datetime.strftime('%Y-%m-%dT%H:%M:%S') for t in times]

            table.add_column(Column(times, name='discovery_date'))

    return(table)

# Provides table names as keyword arguments
def get_table_names(event):

    base_path = pathlib.Path(__file__).parent.resolve()
    data_dir = os.path.join(base_path, 'data')
    base_event_dir = os.path.join(data_dir, event)

    if not os.path.exists(base_event_dir):
        os.makedirs(base_event_dir)

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

    # Make all table names full file names
    for key in tbl.keys():
        tbl[key] = os.path.join(base_event_dir, tbl[key])

    tbl['event_dir']=base_event_dir
    tbl['data_dir']=data_dir
    tbl['base_dir']=base_path

    return(tbl)

def parse_kwargs(args):

    event = args.gw_event
    redo = args.redo
    map_type = args.map_type

    source_table = get_candidate_table(args.sources)
    constraints = add_constraints(args)

    kwargs = {'event': event,
                'meta': [],
                'reference': None,
                'reference_name': '',
                'yse': source_table,
                'redo': redo,
                # Filters to avoid for photometric analysis of GW candidates
                'avoid_filters': ['G','B','V','U','UVM2','UVW2',
                                  'UVW1','orange','cyan','R'],
                'redshift': {
                    'methods': args.redshift_methods,
                    'use_name': 'use_redshift',
                },
                # astcheck observatory code - see: 
                # https://minorplanetcenter.net//iau/lists/ObsCodes.html
                # Default is 304=Las Campanas Observatory
                'obscode': args.obscode,
                'tns': {'api_key': args.tns_api_key,
                        'bot_name': args.tns_bot_name,
                        'bot_id': args.tns_bot_id,},
                # Table name to pull for photometry on YSE-PZ
                'photyse': args.yse_phot,
                'shibboleth': args.shibboleth,
                # Analysis procs that don't require crossmatching
                'no_crossmatch': ['import candidates', 'time', 'prob',
                    'redshift','class','photometry'],
                'outtable': {
                    'ra_precision': 3,
                    'dec_precision': 2,
                    'sort': ['name_len','name'],
        }
    }

    # Parsing GW event data
    if event is not None:
        tbl = get_table_names(event)
        url=f'https://gracedb.ligo.org/api/superevents/{event}/files/{map_type}.fits.gz'
        event_file = os.path.join(tbl['event_dir'], event+f'_{map_type}.pkl')
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

        kwargs['alert_time']=Time(header['DATE-OBS'])
        kwargs['gw_map_url']=url
        kwargs['gw_map']=prob
        kwargs['gw_map_table']=table
        kwargs['gw_map_header']=header
    else:
        event = 'default'
        tbl = get_table_names(event)

    kwargs['event']=event
    kwargs.update(tbl)
    kwargs.update(constraints)

    return(kwargs)


def parse_steps(args):

    import utilities.utilities as util

    out_steps = []
    for step in args.steps:

        step = step.lower()

        if 'import' in step:
            out_steps.append(util.import_candidates)
        elif 'time' in step:
            out_steps.append(util.check_time)
        elif 'prob' in step:
            out_steps.append(util.check_prob)
        elif 'mpc' in step:
            out_steps.append(util.check_mpc)
        elif 'astcheck' in step:
            astcheck_dependency_check()
            out_steps.append(util.check_astcheck)
        elif 'gaia' in step:
            out_steps.append(util.check_gaia)
        elif 'class' in step:
            out_steps.append(util.check_class)
        elif 'ps1dr2' in step:
            out_steps.append(util.check_ps1dr2)
        elif 'asassn' in step:
            out_steps.append(util.check_asassn)
        elif 'redshift' in step:
            out_steps.append(util.check_redshift)
        elif 'photometry' in step:
            out_steps.append(util.check_photometry)
        elif 'tns' in step:
            out_steps.append(util.check_tns)
        elif 'yse' in step:
            out_steps.append(util.check_yse)
        else:
            raise Exception(f'ERROR: unrecognized step {step}!')

    return(out_steps)

def add_options():
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("sources", type=str,
        help="Source of transient/candidate data to analyze.  Can be a path to file or a URL.")
    parser.add_argument("--gw-event", type=str, default=None,
        help="Input event name.  Must be parseable from GraceDB "+\
        "(required for probability cut).")
    parser.add_argument("--cut-time", type=str, default='',
        help="Override event time with input value (usually parsed from Grace DB).")
    parser.add_argument("-r", "--redo", default=False, action='store_true',
        help="Redo all steps if they have not yet been run.")
    parser.add_argument("-s", "--shibboleth", type=str,
        default='/home/ckilpatrick/scripts/shibboleth',
        help="Shibboleth file for sub-routines that require authentication.")
    parser.add_argument("-m", "--map-type", type=str,
        default="bayestar",
        help="Map type for analysis (bayestar should always exist).")
    # Specific steps to be run by candidate analysis pipeline
    parser.add_argument("--steps", nargs='+', type=str,
        default=["import","time","prob","mpc","gaia",
        "class","ps1dr2","asassn","redshift"], 
        help="List of steps to check before validating a candidate.")
    # Method-specific parameters
    parser.add_argument("--max-date", default=5, type=float,
        help="Maximum days from GW discovery to analyze candidates (in days).")
    parser.add_argument("--max-prob", default=0.95, type=float,
        help="Maximum cumulative probability (from 0-1) to consider for valid candidates.")
    parser.add_argument("--distance-distribution", default=0.95, type=float,
        help="Distance distribution (in sigma; 0-infinity) to consider for candidate analysis.")
    parser.add_argument("-z", "--redshift", nargs=2, type=float, default=[-1., 0.35],
        help="Redshift range in which to search for candidates.")
    # Apparent radius cuts for catalog crossmatching
    parser.add_argument("--mpc-radius", default=20.0, type=float,
        help="Radius for minor planet checker (in arcsec).")
    parser.add_argument("--astcheck-radius", default=20.0, type=float,
        help="Radius for astcheck minor planet code (in arcsec).")
    parser.add_argument("--gaia-radius", default=2.0, type=float,
        help="Radius for crossmatching to Gaia stars (in arcsec).")
    parser.add_argument("--asassn-radius", default=2.0, type=float,
        help="Radius for crossmatching to ASASSN variables (in arcsec).")
    parser.add_argument("--ps1-radius", default=1.0, type=float,
        help="Radius for crossmatching to PS1 DR2 stars (in arcsec).")
    parser.add_argument("--tns-radius", default=3.0, type=float,
        help="Radius for crossmatching to TNS transients (in arcsec).")
    parser.add_argument("--yse-radius", default=3.0, type=float,
        help="Radius for crossmatching to YSE transients (in arcsec).")
    parser.add_argument("--gal-radius", default=60.0, type=float,
        help="Radius for crossmatching to candidate host galaxies (in arcsec).")
    parser.add_argument("--max-gal-radius","--max-galaxy-radius", 
        default=3.0e5, type=float,
        help="Maximum separation from a galaxy with a known distance/redshift"+\
        " (to rule out chance coincidence; in parsec).")
    parser.add_argument("--phot-bright", default=None, type=float,
        help="For candidates with known distance, maximum apparent magnitude.")
    parser.add_argument("--phot-min-decline", default=None, type=float,
        help="Minimum decline rate (mag/day) of candidates to be considered.")
    parser.add_argument("--phot-color", default=None, type=float,
        help="Minimum color (g-r mag) of candidates to be considered.")
    # Specific parameters for search and redshift methods
    parser.add_argument("--yse-phot", type=int, default=None,
        help="ID of YSE-PZ photometry SQL query.  Default is not to use one.")
    parser.add_argument("--redshift-methods", nargs='+',
        default=["spec","ps1strm","2mpz","legacy"],
        help="List of sources to check for host redshifts [spec,ps1strm,2mpz,legacy]")
    parser.add_argument("--tns-api-key", type=str, default=None,
        help="TNS API key - required for TNS queries.")
    parser.add_argument("--tns-bot-name", type=str, default="YSE_progenitors",
        help="TNS API key - required for TNS queries.")
    parser.add_argument("--tns-bot-id", type=str, default="97993",
        help="TNS API key - required for TNS queries.")
    parser.add_argument("--latex", default=False, action='store_true',
        help="Output latex-formatted candidate file.")
    parser.add_argument("--candidate-format", default=False,action='store_true',
        help="Output candidate-formatted file.")
    parser.add_argument("--obscode", default="304", type=str,
        help="Observatory code for astcheck functionality.  For a list, see "+\
        "https://minorplanetcenter.net//iau/lists/ObsCodes.html.")

    args = parser.parse_args()

    return(args)
