import warnings
warnings.filterwarnings('ignore')
from astropy.io import ascii
from astropy.table import Table, Column, vstack, unique, hstack
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck15 as cosmo
from astropy import units as u
from bs4 import BeautifulSoup
from scipy.optimize import curve_fit
import copy, sys, requests, os, time, pandas, progressbar, pickle, glob
import healpy as hp
import numpy as np
import astropy_healpix as ah
from ligo.skymap import distance as lvc_distance
from astroquery.ned import Ned
from astroquery.vizier import Vizier

Vizier.ROW_LIMIT = -1

from contextlib import contextmanager
@contextmanager
def suppress_stdout():
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stdout

with suppress_stdout():
    from astropy.time import Time
    from mwdust.util import extCurves
    from dustmaps.sfd import SFDQuery
    import dustmaps

try:
    sfd = SFDQuery()
except:
    dustmaps.sfd.fetch()
    sfd = SFDQuery()

def lin_decline(x, a, b):
    return(a+b*x)

def is_number(val):
    try:
        val = float(val)
        if not np.isnan(val):
            return(True)
    except:
        pass

    return(False)

def sanitize_table(table, reformat_columns=False, reformat_z=False,
    remove_columns=False):

    remove=[]
    if len(table)>0:
        for key in table.keys():
            var = table[key][0]
            if isinstance(var, (list, np.ndarray)):
                remove.append(key)

    for r in remove:
        table.remove_column(r)

    for key in table.keys():
        # Avoid deleting keys that are already lower case
        if key==key.lower().strip():
            continue

        if '_absmag' in key or '_var' in key:
            continue

        if key.lower().strip() not in table.keys():
            table.rename_column(key, key.lower().strip())
        else:
            table.remove_column(key)

    for key in ['masterid']:
        if key in table.keys() and 'objid' not in table.keys():
            table.rename_column(key, 'objid')

    if 'cz' in table.keys() and 'z' not in table.keys():
        # Mask table to make sure values are all numbers
        mask = [is_number(val) for val in table['cz']]
        table = table[mask]
        col = [float(val)/2.998e5 for val in table['cz'].data]
        table.add_column(Column(col, name='z'))

    if 'e_cz' in table.keys() and 'z_err' not in table.keys():
        mask = [is_number(val) for val in table['e_cz']]
        table = table[mask]
        col = [float(val)/2.998e5 for val in table['e_cz'].data]
        table.add_column(Column(col, name='z_err'))

    # Fix issue that sometimes occurs with "name" column
    for key in table.keys():
        if 'name' in key:
            table.rename_column(key, 'name')


    # Try to put the table into a standard format with a name, ra, dec, time,
    # spectroscopic class, and redshift column
    if reformat_columns:

        key_mapping = {
            'ra': ['transient_ra','ramean','raj2000', 'right ascension'],
            'dec': ['transient_dec','decmean','dej2000', 'declination'],
            'classification': ['spec_class'],
            'discovery_date': ['disc_date', 'date'],
            'name': ['objid','objectid','uniquePspsP2id','id'],
            'z': ['z_phot_median','z_phot'],
            'z_err': ['z_phot_std','z_photerr']
        }

        if reformat_z:
            del key_mapping['z']
            del key_mapping['z_err']
            key_mapping['redshift']=['z']

        for key in key_mapping.keys():
            if key in table.keys():
                continue
            else:
                for other_key in key_mapping[key]:
                    if other_key in table.keys():
                        table.rename_column(other_key, key)
                        break

        if remove_columns:
            for key in table.keys():
                if key not in key_mapping.keys():
                    table.remove_column(key)

    # Sanitize masks
    for key in table.keys():
        if '_mask' in key:
            newdata = []
            for val in table[key].data:
                if 'false' in str(val).lower():
                    newdata.append(False)
                else:
                    newdata.append(True)

            table.remove_column(key)
            table.add_column(Column(newdata, name=key))

    return(table)

def get_legacy_file(coord, *args):
    ra = coord.ra.degree
    dec = coord.dec.degree
    ra_min = int(np.floor(ra/10)) * 10
    ra_max = int(np.ceil(ra/10)) * 10
    sign = 'p'
    if dec < 0: sign='m'
    dec_min = np.abs(int(np.floor(dec/5))) * 5
    dec_max = np.abs(int(np.ceil(dec/5))) * 5
    file_fmt = 'sweep-{ra_min}{sign}{dec_min}-{ra_max}{sign}{dec_max}.fits'
    file = file_fmt.format(ra_min=str(ra_min).zfill(3),
        ra_max=str(ra_max).zfill(3), dec_min=str(dec_min).zfill(3),
        dec_max=str(dec_max).zfill(3), sign=sign)
    return(file)

def search_legacy_files(coord, *args, radius=30 * u.arcsec,
    base_dir='/data/Legacy/photoz'):
    ra = coord.ra.degree
    dec = coord.dec.degree
    rad = radius.to_value('arcsec') / 3600.
    ra_lim = np.array([-rad, rad])*np.cos(dec * np.pi/180.) + ra
    dec_lim = np.array([-rad, rad]) + dec
    files = []
    outtables = []
    for r in ra_lim:
        rsearch = r
        if r < 0:
            r_limits = [r+360.0, 360.0]
            rsearch = r + 360.0
        elif r > 360:
            r_limits = [0., r-360.0]
            rsearch = r - 360.0
        else:
            r_limits = ra_lim
        for d in dec_lim:
            dsearch = d
            if d < -90.0:
                d_limits = [-90.0, dec_lim[1]]
                dsearch = -89.99
            elif d > 90.0:
                d_limits = [dec_lim[0], 90.]
                dsearch = 89.99
            else:
                d_limits = dec_lim
            file = get_legacy_file(SkyCoord(rsearch, dsearch, unit='deg'))

            fullfile = os.path.join(base_dir, file)
            fullfile_pz = fullfile.replace('.fits','-pz.fits')

            if (file not in files and os.path.exists(fullfile) and
                os.path.exists(fullfile_pz)):
                files.append(file)

                table = Table.read(fullfile, hdu=1)
                mask = ((table['RA'] < r_limits[1]) &
                    (table['RA'] > r_limits[0]) &
                    (table['DEC'] > d_limits[0]) &
                    (table['DEC'] < d_limits[1]))

                pztable = Table.read(fullfile_pz, hdu=1)
                outtables.append(hstack([table[mask], pztable[mask]]))

    if len(outtables)>0:
        outtable = vstack(outtables)

        for key in outtable.keys():
            outtable.rename_column(key, key.lower())
        return(outtable)

    else:
        return(None)

def get_ps1strm(coord, *args, radius=300 * u.arcsec, **kwargs):

    ra = coord.ra.degree ; dec = coord.dec.degree
    if dec < -30.0:
        return(None)

    ra = '-'+str(coord.ra.degree)
    if dec < 0:
        dec = '+'+str(np.abs(coord.dec.degree))
    else:
        dec = '-'+str(coord.dec.degree)

    catalog_dir = '/data/PS1STRM/photoz'
    catalogs = sorted(glob.glob(os.path.join(catalog_dir, '*strm*')))

    limits = []
    for catalog in catalogs:
        dec_str = catalog.split('_')[4]
        demin, demax = dec_str.split('-')

        if demin.startswith('m'):
            demin = -1.0 * float(demin[1:])
        else:
            demin = float(demin[1:])

        if demax.startswith('m'):
            demax = -1.0 * float(demax[1:])
        else:
            demax = float(demax[1:])

        limits.append((float(demin), float(demax)))

    idx=0
    for lim in limits:
        if (float(dec) > float(lim[0])) and (float(dec) < float(lim[1])):
            break
        else:
            idx += 1

    if idx>=len(catalogs):
        return(None)

    catalog = catalogs[idx]

    file = 'tmp'

    cmd = 'zcat {0} | '.format(catalog)
    cmd += 'awk -F\',\' \'{if (($3 '
    cmd += '{0} )^2+($4 '.format(ra)
    cmd += '{0} )^2 <'.format(dec)
    cmd += ' ( {0} )^2)'.format(radius.to_value('degree'))
    cmd += ' print}\''
    cmd += ' > {0}'.format(file)

    os.system(cmd)

    catnames = 'objID,uniquePspsOBid,raMean,decMean,l,b,class,prob_Galaxy,'+\
    'prob_Star,prob_QSO,extrapolation_Class,cellDistance_Class,cellID_Class,'+\
    'z_phot,z_photErr,z_phot0,extrapolation_Photoz,cellDistance_Photoz,'+\
    'cellID_Photoz'
    names = catnames.split(',')

    try:
        table = ascii.read(file, names=names)
        if os.path.exists('tmp'): os.remove('tmp')
        return(table)
    except:
        if os.path.exists('tmp'): os.remove('tmp')
        return(None)

def get_css(ra, dec, radius):

        url = 'http://nunuku.caltech.edu/cgi-bin/getcssconedb_release_img.cgi'
        form_data={'RA': ra,
                   'Dec': dec,
                   'Rad': radius,
                   'DB': 'photcat',
                   'SHORT': 'short',
                   'OUT': 'csv'}

        r = requests.post(url, data=form_data)

        regex = r'result\_web\_file([a-zA-Z0-9]+)\.csv'
        match = re.search(regex, r.text)

        if match:
            filename = match.group(0)
            file_url = 'http://nunuku.caltech.edu/DataRelease/upload/'+filename
            table = Table.from_html(file_url)

            with suppress_stdout:
                try:
                    table = Table.read(file_url)
                    table = sanitize_table(table)
                    return(table)
                except:
                    return(None)
        else:
            return(None)

def check_ps1dr2(coord, *args, **kwargs):

    radius = kwargs['search_radius']['ps1dr2'].to_value('degree')/36.
    region = '{0} {1}'.format(coord.ra.degree, coord.dec.degree)

    from astroquery.mast import Catalogs
    catalog_data = Catalogs.query_region(region, radius=radius,
        catalog='Panstarrs',data_release='dr2',table='detections')

    catalog_data = sanitize_table(catalog_data, reformat_columns=True)

    return(catalog_data)

def check_mpc(coord, transient_time, *args, **kwargs):
        table_output = []
        day = transient_time.strftime('%d')
        dayfloat = transient_time.mjd % 1.0
        day = day + '.' + str(dayfloat)[2:]
        # Hack AF
        url = 'https://minorplanetcenter.net/cgi-bin/mpcheck.cgi'
        form_data = {"year":"%s" % transient_time.strftime("%Y"),
            "month":"%s" % transient_time.strftime("%m"),
            "day":"%s" % day,
            "which":"pos",
            "ra":"%s" %
                coord.ra.to_string(unit="hour",pad=True,decimal=False,sep=" "),
            "decl":"%s" %
                coord.dec.to_string(unit="deg",pad=True,decimal=False,sep=" "),
            "TextArea":"",
            "radius":"%s" % 30,
            "limit":"%s" % 24,
            "oc":"%s" % 500,
            "sort":"d",
            "mot":"h",
            "tmot":"s",
            "pdes":"u",
            "needed":"f",
            "ps":"n",
            "type":"p"}

        r = requests.post(url, data=form_data)
        soup = BeautifulSoup(r.text,'lxml')
        pre = soup.find("pre")
        if pre is None:
            return(None)
        else:
            data = []
            for row in [row for row in pre.contents[-1].split("\n")[3:-1]]:
                coord = SkyCoord(row[25:36], row[36:47], unit=(u.hour, u.deg))
                data.append([row[9:25],coord.ra.degree,coord.dec.degree])
            table = Table(list(map(list, zip(*data))),names=('name','ra','dec'))
            return(table)

def check_gaia(coord, check_time, *args, **kwargs):
    radius = u.Quantity(30.0, u.arcsec)
    with suppress_stdout():
        from astroquery.gaia import Gaia
        job = Gaia.cone_search_async(coord, radius, verbose=False)
    r = job.get_results()

    # Do cut on parallax.  Must be >3 sigma
    match = r['parallax']/r['parallax_error'] > 3.0
    if len(r[match])>0:
        return(r[match])
    else:
        return(None)

def check_asassn(coord, transient_time, *args, **kwargs):
    # Max radius is 10.0
    radius = 10.0
    url='https://asas-sn.osu.edu/photometry?utf8=%E2%9C%93&ra={ra}&dec={dec}'+\
        '&radius={rad}&vmag_min=&vmag_max=&epochs_min=&epochs_max=&rms_min=&'+\
        'rms_max=&sort_by=raj2000'
    url = url.format(ra=coord.ra.degree, dec=coord.dec.degree, rad=radius)
    r = requests.get(url)
    soup = BeautifulSoup(r.text)
    tab = soup.find('table')
    try:
        df = pandas.read_html(str(tab))
    except ValueError:
        return(None)
    table = Table.from_pandas(df[0])

    table = sanitize_table(table, reformat_columns=True)

    return(table)

def search_spectroscopic_redshift(coord, *args, radius=300*u.arcsec, **kwargs):

    result_table = Ned.query_region(coord, radius=radius)

    mask = (result_table['Redshift Flag']=='SPEC') & (result_table['Type']=='G')
    result_table = result_table[mask]

    result_table = result_table['Object Name','RA','DEC','Redshift']

    result_table.rename_column('RA', 'ra')
    result_table.rename_column('DEC', 'dec')
    result_table.rename_column('Redshift', 'z')
    result_table.rename_column('Object Name', 'name')

    coords = np.array([parse_coord(r['ra'], r['dec']) for r in result_table])
    sep = np.array([coord.separation(c).arcmin for c in coords])

    result_table.add_column(Column(sep, name='separation'))

    z_err = result_table['z'].data*1.0e-3
    result_table.add_column(Column(z_err, name='z_err'))

    return(result_table)

def get_2mass_redshift(coord, *args, radius=300*u.arcsec, **kwargs):

    v = Vizier(columns=['ID','RAJ2000', 'DEJ2000','cz', 'e_cz'])
    result_table = v.query_region(coord, radius=radius,
        catalog='J/ApJS/199/26/table3')

    if len(result_table)==0:
        return(None)
    else:
        result_table = result_table[0]

    z = result_table['cz'].data/(2.998e5)
    e_z = result_table['e_cz'].data/(2.998e5)

    result_table.add_column(Column(z, name='z'))
    result_table.add_column(Column(e_z, name='z_err'))

    result_table.rename_column('RAJ2000','ra')
    result_table.rename_column('DEJ2000','dec')
    result_table.rename_column('ID','name')

    coords = SkyCoord(result_table['ra'], result_table['dec'], unit='deg')

    sep = coord.separation(coords)

    result_table.add_column(Column(sep, name='separation'))

    return(result_table)


def build_meta_table(table, typ, r=30.0 * u.arcsec):
    methods={'legacy': search_legacy_files, 'ps1strm': get_ps1strm,
        'spec': search_spectroscopic_redshift, '2mpz': get_2mass_redshift}
    meta_table = None
    bar = progressbar.ProgressBar(maxval=len(table))
    bar.start()
    for i,row in enumerate(table):
        bar.update(i+1)
        time.sleep(0.1)

        coord = parse_coord(row['ra'], row['dec'])

        if not meta_table:
            meta_table = methods[typ](coord, radius=r)
        else:
            subtable = methods[typ](coord, radius=r)

            if subtable and len(subtable)>0:
                meta_table = vstack([subtable, meta_table])
                for key in meta_table.keys():
                    try:
                        meta_table = meta_table[~meta_table[key].mask]
                    except AttributeError:
                        pass

                meta_table = unique(meta_table)

    bar.finish()

    return(meta_table)

def add_name_len(table):
    if 'name' not in table.keys():
        return(table)

    try:
        name_len = [len(row['name']) for row in table]
        table.add_column(Column(name_len, name='name_len'))
        return(table)
    except:
        pass

    return(table)

def format_redshift(row, **kwargs):
    out = '--'
    use_name = ''
    if kwargs['redshift']['use_name'] in row.colnames:
        use_name = kwargs['redshift']['use_name']
    else:
        return(out)

    if 'mpc_mask' in row.colnames and row['mpc_mask']:
        return('--')

    if row[use_name]!='--':
        use = row[use_name]
        z_fmt = '{0}$\\pm${1} ({2})'

        z = '' ; z_err = ''
        if use+'_z' in row.colnames and use+'_z_err' in row.colnames:
            z = '%7.5f'%float(row[use+'_z'])
            z_err = '%7.5f'%float(row[use+'_z_err'])
        else:
            return(out)

        source=''
        if use=='spec': source='s'
        if use=='ps1strm': source='PS1'
        if use=='legacy': source='LDR10'
        if use=='2mpz': source='2MRS'

        if not source: return(out)

        out = z_fmt.format(z, z_err, source).strip()

        return(out)
    else:
        return(out)

def format_coord(row, typ, **kwargs):
    coord = parse_coord(row['ra'], row['dec'])
    if coord:
        p = kwargs['outtable'][typ+'_precision']
        data = coord.to_string(style='hmsdms', precision=p, sep=':').split()
        if typ=='ra': return(data[0])
        if typ=='dec': return(data[1])

    return('--')

def format_variability(row, filt, **kwargs):
    delta_mag = '(PHOT; $\\Delta_{0}$={1}$\\pm${2})'
    delta = '%2.2f'%row[filt+'_var']
    err = '%2.2f'%row[filt+'_var_err']
    dat = delta_mag.format('{'+filt+'}', delta, err)

    return(dat)

def format_note(row, unique_filters, **kwargs):
    note = ''
    if 'class_mask' in row.colnames and row['class_mask']:
        note='(SN) '+row['classification'].replace('SN','').strip()
    elif 'mpc_mask' in row.colnames and row['mpc_mask']:
        note='(MP) '+row['mpc']
        z='--'
    elif 'redshift_mask' in row.colnames and row['redshift_mask']:
        note='(Z)'
    elif 'absmag_mask' in row.colnames and row['absmag_mask']:
        note = '(PHOT; bright)'
    elif 'var_mask' in row.colnames and row['var_mask']:
        use_filt = ''
        for filt in unique_filters:
            if row[filt+'_var_mask']: use_filt = filt
        if use_filt:
            note = format_variability(row, use_filt, **kwargs)

    note = note.strip()

    return(note)

def format_var(row, unique_filts, **kwargs):
    var='--'

    if row['mpc_mask']:
        return(var)

    # Pick which filter to use for absolute magnitude
    use_filt = '--'
    for filt in unique_filts:
        if row[filt+'_absmag_mask']:
            use_filt = filt

    if use_filt=='--':
        data = [[float(row[filt+'_absmag_date']),
                 float(row[filt+'_absmag']),
                 float(row[filt+'_absmag_err'])] for filt in unique_filts]
        use_data = []
        for i,val in enumerate(data):
            if all([(not np.isnan(v) and v is not None) for v in val]):
                val.append(i)
                use_data.append(val)
        if use_data:
            use_data = sorted(use_data, key=lambda v: v[0])
            use_filt = unique_filts[use_data[0][3]]

    if use_filt=='--' and row['absmag_filt']!='--':
        use_filt = row['absmag_filt']

    filt = use_filt
    if not filt=='--':
        if (filt+'_absmag' in row.colnames and filt+'_absmag_err' and
            filt+'_absmag_date' in row.colnames):

            if all([not np.isnan(row[filt+'_absmag'+typ])
                for typ in ['','_err','_date']]):

                var_fmt = '${0}$=$-${1}$\\pm${2} ({3} d)'
                mag = '%2.2f'%row[filt+'_absmag']
                magerr = '%2.2f'%row[filt+'_absmag_err']
                date = '%2.1f'%float(row[filt+'_absmag_date'])
                var = var_fmt.format(filt, mag, magerr, date)

    return(var)

def output_latex_table(table, **kwargs):
    varnams = ['name','ra','dec','prob','date','z','var','note']
    table = add_name_len(table)
    ufilts = table.meta['unique_filters']
    if all([n in table.keys() for n in kwargs['outtable']['sort']]):
        table.sort(kwargs['outtable']['sort'])

    # Mask the table to only include things that pass time and prob cuts
    if 'prob_mask' in table.keys() and 'time_mask' in table.keys():
        mask = ~table['prob_mask'] & ~table['time_mask']
        table = table[mask]

    all_data = []
    for row in table:

        # Format coordinate
        data = {}
        for var in varnams:
            var = var.strip()
            if var=='name': data[var]=row[var]
            if var=='ra' or var=='dec': data[var]=format_coord(row,var,**kwargs)
            if var=='date': data[var] = '%7.5f'%Time(row['discovery_date']).mjd
            if var=='prob': data[var] = '%1.4f' % row['2d_probability']
            if var=='z': data[var] = format_redshift(row, **kwargs)
            if var=='var': data[var] = format_var(row, ufilts, **kwargs)
            if var=='note': data[var] = format_note(row, ufilts, **kwargs)

        all_data.append(data)

    lens = [10]*len(varnams)
    for i,var in enumerate(varnams):
        var = var.strip()
        max_length = np.max([len(str(dat[var])) for dat in all_data])
        lens[i] = max_length

    out_fmt = ' & '.join(['{'+v+': <'+str(l)+'}' for v,l in zip(varnams, lens)])

    hdr={}
    for key in varnams: hdr[key]=key.upper()
    print('\n\n'+out_fmt.format(**hdr)+'\n\n')
    for data in all_data:
        output = out_fmt.format(**data)
        print(output + '\\\\')

def get_dist_sigma(dist_map, sigma_map, pix):
    d = dist_map[pix]
    s = sigma_map[pix]
    mean, stddev, norm = lvc_distance.parameters_to_moments(d, s)
    return(mean, stddev)

# Downloads an essentially complete list (look at SQL query params) of all
# supernovae in the YSE PZ data base
def download_yse(yse):

    url = 'https://ziggy.ucolick.org/yse/explorer/{0}/download?format=csv'
    url = url.format(yse)

    data = requests.get(url)
    table = ascii.read(data.text)

    for key in table.keys():
        table.rename_column(key, key.lower())

    use_key = None
    for key in ['spec_class','classification','class']:
        if key in table.keys():
            use_key = key
            break

    if use_key is not None:
        mask = table[use_key]!='FRB'
        table = table[mask]

    for key in table.keys():
        if 'transient_ra' in key.lower(): table.rename_column(key, 'ra')
        if 'transient_dec' in key.lower(): table.rename_column(key, 'dec')
        if 'disc_date' in key.lower(): table.rename_column(key, 'discovery_date')
        if 'spec_class' in key.lower(): table.rename_column(key, 'classification')
        if 'name' in key.lower(): table.rename_column(key, 'name')

    return(table)

def parse_coord(ra, dec):
    def check_coord(ra, dec):
        if ':' in str(ra) and ':' in str(dec):
            coord = SkyCoord(ra, dec, unit=(u.hour, u.deg))
            return(coord)
        else:
            try:
                ra = float(ra) ; dec = float(dec)
                coord = SkyCoord(ra, dec, unit=(u.deg, u.deg))
                return(coord)
            except:
                return(None)

    if (isinstance(ra, (list, np.ndarray, Column)) and
        isinstance(dec, (list, np.ndarray, Column))):
        coords = np.array([check_coord(r,d) for r, d in zip(ra, dec)])
        return(coords)
    else:
        return(check_coord(ra, dec))

def import_candidates(table, **kwargs):

    table_name = kwargs['candidate']

    if not os.path.exists(table_name) or kwargs['redo']:
        table = download_yse(kwargs['yse'])
        table = sanitize_table(table, reformat_columns=True, reformat_z=True)
        table = table['name','ra','dec','discovery_date']

        # If we want to import additional candidates from a separate file
        if kwargs['internal'] and os.path.exists(kwargs['internal']):
            additional_table = ascii.read(kwargs['internal'])
            additional_table = sanitize_table(additional_table)

            for row in additional_table:
                coord = parse_coord(row['ra'], row['dec'])

                coords = [parse_coord(r['ra'], r['dec']) for r in table]
                coords = np.array(coords)

                if coord:
                    radius = kwargs['search_radius']['asassn']
                    match = np.array([coord.separation(c) < radius for c in coords])

                    if len(table[match])>0:
                        continue
                    else:
                        t = Time(row['mjd'], format='mjd')
                        table.add_row([row['name'], coord.ra.degree,
                            coord.dec.degree,' '.join(t.fits.split('T'))])

        if kwargs['gw_map_url']:
            if kwargs['gw_map'] is None:
                filename = download_file(kwargs['gw_map_url'], cache=True)
                prob, header = hp.read_map(filename, h=True)
                gw_table = Table.read(kwargs['gw_map_url'])
            else:
                gw_map = kwargs['gw_map']
                header = kwargs['gw_map_header']
                gw_table = kwargs['gw_map_table']

            probability=[]
            nside = header['NSIDE']

            sorted_pixels = np.flip(np.argsort(gw_map))
            sorted_prob = np.flip(sorted(gw_map))

            def ci(level, sorted_prob):
                csum = 0
                c = 0
                index = 0
                while csum < level:
                    csum += sorted_prob[index]
                    c = sorted_prob[index]
                    index += 1
                return csum, c, index

            c99sum, c99, index = ci(0.99, sorted_prob)

            print('Largest probability pixel',gw_map[sorted_pixels[0]])
            print('Smallest probability pixel',gw_map[sorted_pixels[-1]])

            area_per_pix = hp.pixelfunc.nside2pixarea(header['NSIDE'],
                degrees=True)

            total_prob = 0.0
            total_area = 0.0
            least_prob = 0.0
            for i,pix in enumerate(sorted_pixels):
                if total_prob > kwargs['probability']:
                    least_prob = gw_map[pix]
                    break
                else:
                    total_prob += gw_map[pix]
                    total_area += area_per_pix

            total_prob_str = '%1.6f'%total_prob
            total_area_str = '%5.4f'%total_area
            print(f'Total area for {total_prob_str} probability '+\
                f'is {total_area_str} deg^2')

            # Add this information to the table metadata
            table.meta['gw_probability']=kwargs['probability']
            table.meta['gw_probability_area']=total_area

            # Create a list of cumulative sum of the value of pixels in hp map
            cumulative = np.zeros(len(sorted_pixels))
            approx = 0.999999999999
            for i in np.arange(len(sorted_pixels)):
                # Don't need to go beyond very low prob
                if cumulative[sorted_pixels[i-1]] > approx:
                    cumulative[sorted_pixels[i:]] = 1.0
                    break
                if i==0:
                    cumulative[sorted_pixels[i]] = gw_map[sorted_pixels[i]]
                else:
                    idx0 = sorted_pixels[i-1]
                    idx1 = sorted_pixels[i]
                    cumulative[idx1] = cumulative[idx0] + gw_map[idx1]

            print(f'It took {i} steps to calculate up to {approx} probability')

            probability = []
            for row in table:
                coord = parse_coord(row['ra'], row['dec'])

                ra = coord.ra.degree
                dec = coord.dec.degree

                theta = 0.5 * np.pi - np.deg2rad(dec)
                phi = np.deg2rad(ra)

                ipix = hp.ang2pix(header['NSIDE'], theta, phi)

                probability.append(cumulative[ipix])

            table.add_column(Column(probability, name='2d_probability'))

            table = add_distance_data(table, gw_map=gw_table, header=header,
                table_name=kwargs['candidate'])

        # Write out the table with basic information
        table.write(table_name, format='ascii', overwrite=True)

    else:
        table = ascii.read(table_name)

    return(table)

def check_class(table, **kwargs):

    # Add extra classifications to table
    classtable = download_yse('45')

    classification = []
    for row in table:
        match = classtable['name']==row['name']
        if len(classtable[match])>0:
            classification.append(classtable[match][0]['classification'])
        else:
            classification.append('--')

    if 'classification' in table.keys():
        table.remove_column('classification')

    table.add_column(Column(classification, name='classification'))

    # Reset classification mask column if it does not exist
    if 'class_mask' in table.keys():
        table.remove_column('class_mask')

    class_col = Column([any([c in str(row['classification'])
        for c in ['Ia','II','Ib','Ic','TDE','AGN']])
        for row in table], name='class_mask')
    table.add_column(class_col)

    return(table)

def add_distance_data(table, gw_map='', header={}, table_name=''):

    if ('DL' not in table.keys() or 'DL_ERR' not in table.keys() and gw_map):

        dist_map = gw_map['DISTMU']
        sigma_map = gw_map['DISTSIGMA']

        nside = header['NSIDE']

        area_per_pix = hp.pixelfunc.nside2pixarea(nside, degrees=True)

        # Now add candidate cumulative probability to the table
        distance = [] ; dist_sigma = []
        for row in table:
            coord = parse_coord(row['ra'], row['dec'])

            ra = coord.ra.degree
            dec = coord.dec.degree

            theta = 0.5 * np.pi - np.deg2rad(dec)
            phi = np.deg2rad(ra)

            ipix = hp.ang2pix(header['NSIDE'], theta, phi)

            dist, sigma = get_dist_sigma(dist_map, sigma_map, ipix)

            distance.append(dist)
            dist_sigma.append(sigma)

        table.add_column(Column(distance, name='DL'))
        table.add_column(Column(dist_sigma, name='DL_ERR'))

        if table_name:
            table.write(table_name, overwrite=True, format='ascii')

    return(table)

def redshift_comparison(table, key, z_sigma_factor=1.0, dl_sigma_factor=2.58):

    if key+'_z' not in table.keys() or key+'_z_err' not in table.keys():
        return(table)

    mask=[]
    for row in table:
        try:
            # Assuming key is a redshift
            z = float(row[key+'_z'])
            z_err = float(row[key+'_z_err'])
            dl = float(row['dl'])
            dl_err = float(row['dl_err'])
        except ValueError:
            mask.append(False)
            continue

        z_range=(z-z_sigma_factor*z_err, z+z_sigma_factor*z_err)
        distance_range = (cosmo.luminosity_distance(z_range[0]).value,
            cosmo.luminosity_distance(z_range[1]).value)

        dl_range = (dl-dl_sigma_factor*dl_err, dl+dl_sigma_factor*dl_err)

        if distance_range[0] > dl_range[1] or distance_range[1] < dl_range[0]:
            mask.append(True)
        else:
            mask.append(False)

    if key+'_mask' in table.keys():
        table.remove_column(key+'_mask')

    table.add_column(Column(mask, name=key+'_mask'))
    table = sanitize_table(table)

    return(table)

def host_crossmatch(row, reference, radius=300.0 * u.arcsec,
    limit=3.0e5 * u.parsec):

    if len(reference)==0:
        return(None)

    check_keys = ['ra', 'dec', 'z', 'z_err']

    for key in check_keys:
        if key not in reference.keys():
            warning = 'WARNING: could not cross reference host galaxies'
            print(warning)
            return(None)

    coord = parse_coord(row['ra'], row['dec'])

    if 'coord' not in reference.keys():
        coords = np.array([parse_coord(r['ra'], r['dec']) for r in reference])
        reference.add_column(Column(coords, name='coord'))

    ang_sep = SkyCoord(reference['coord'].data).separation(coord)
    if 'ang_sep' in reference.keys():
        reference.remove_column('ang_sep')
    reference.add_column(Column(ang_sep, name='ang_sep'))

    match = ang_sep < radius

    if len(reference[match])>0:
        subtable = reference[match]
        unit = subtable['ang_sep'].unit
        sep = [(val['ang_sep'] * unit).to_value('radian') *\
                cosmo.luminosity_distance(val['z']).to_value('parsec')
                for val in subtable]
        separation = Column(sep, name='proj_sep', unit=u.parsec)

        if 'proj_sep' in subtable.keys():
            subtable.remove_column('proj_sep')
        subtable.add_column(separation)

        submatch = subtable['proj_sep'] < limit

        if len(subtable[submatch])>0:
            subtable.sort('proj_sep')
            return(subtable[0])

    return(None)

def format_proc(proc, upper=False):
    try:
        typname = proc.__name__.lower()
        typname = typname.replace('check','')
        typname = typname.replace('get','')
        typname = typname.replace('_',' ')
        typname = typname.strip()
        if upper: typname = typname.upper()
        return(typname)
    except:
        return('')

def add_z(table, ref, procname, add_type, **kwargs):

    #redo = kwargs['redo']
    redo = True

    # We need distance data for the host_crossmatch algorith, so add
    if 'dl' not in table.keys() or 'dl_err' not in table.keys():
        table = add_distance_data(table)

    coltypes = ['ra','dec','z','z_err']

    if any([procname+'_'+typ not in table.keys() for typ in coltypes]) or redo:

        cols = np.array([np.array([np.nan]*len(table))] * len(coltypes))
        bar = progressbar.ProgressBar(max_value=len(table))
        bar.start()
        for i,row in enumerate(table):
            bar.update(i+1)
            match = None
            if add_type=='name': match = ref[ref['name']==row['name']]
            if add_type=='crossmatch': match = host_crossmatch(row, ref)
            for j,typ in enumerate(coltypes):
                if typ and match and len(match)>0:
                    if typ in match.colnames:
                        cols[j,i]=match[typ]
                    else:
                        cols[j,i]=np.nan

        bar.finish()

        for j,typ in enumerate(coltypes):
            if procname+'_'+typ in table.keys():
                table.remove_column(procname+'_'+typ)
            table.add_column(Column(cols[j], name=procname+'_'+typ))

    table = redshift_comparison(table, procname)
    table = sanitize_table(table)

    if kwargs['candidate']:
        table.write(kwargs['candidate'], format='ascii', overwrite=True)

    return(table)

def check_redshift(table, **kwargs):

    for procname in kwargs['redshift']['methods']:
        #if procname in table.keys():
        #    continue

        m='\tREDSHIFT: {proc}'
        print(m.format(proc=procname.upper()))

        method = 'crossmatch'
        #if procname=='spec': method = 'name'

        # Otherwise need to grab reference table and crossmatch
        if os.path.exists(kwargs[procname]):
            reference = ascii.read(kwargs[procname])
            reference = sanitize_table(reference, reformat_columns=True)
        else:
            reference = build_meta_table(table, procname,
                r=kwargs['search_radius']['galaxy'])
            reference = sanitize_table(reference, reformat_columns=True)
            reference.write(kwargs[procname], format='ascii', overwrite=True)

        mask1 = reference['z'] > kwargs['redshift_range'][0]
        mask2 = reference['z'] < kwargs['redshift_range'][1]
        #mask3 = reference['z'] / reference['z_err'] > 1.0

        mask = mask1 & mask2 #& mask3
        reference = reference[mask]

        table = add_z(table, reference, procname, method, **kwargs)

    table = pick_best_redshift(table, **kwargs)

    return(table)

def add_data(table, proc, **kwargs):

    procname = format_proc(proc)
    redo = kwargs['redo']

    # Rapid procs that don't require crossmatching
    if procname in kwargs['no_crossmatch']:
        table = proc(table, **kwargs)
        return(table)

    if procname not in table.keys() or redo:

        reftable = kwargs['reference']
        coltypes = ['match','name','ra','dec']
        coldata = [False, '--', np.nan, np.nan]

        # Load asynchronous data if we need to
        if kwargs['reference_name'] and not reftable:
            if os.path.exists(kwargs['reference_name']):
                rreftable = ascii.read(kwargs['reference_name'])
            else:
                radius = kwargs['search']
                reftable = build_meta_table(table, proc, r=radius)
                reftable.write(kwargs['reference_name'], overwrite=True,
                    format='ascii')

        coldata = []
        bar = progressbar.ProgressBar(max_value=len(table))
        bar.start()

        for i,row in enumerate(table):
            bar.update(i+1)

            discovery_time = Time(row['discovery_date'])
            coord = parse_coord(row['ra'], row['dec'])

            if reftable:
                reference = reftable
            else:
                reference = proc(coord, discovery_time, **kwargs)

            best_match = None
            if reference:
                radius = kwargs['search_radius'][procname]
                coords = parse_coord(reference['ra'], reference['dec'])
                separation = [c.separation(coord).degree for c in coords]
                separation = np.array(separation)
                match =  separation < radius.to_value('degree')

                if len(reference[match])>0:
                    reference.add_column(Column(separation, name='sep'))
                    reference.sort('sep')
                    best_match=reference[0]

            if best_match:
                b=best_match
                name_keys = ['name','DESIGNATION','source_id']
                use_name = ''
                for n in name_keys:
                    if n in b.keys():
                        use_name = n
                if not use_name:
                    print(b)
                    raise Exception('ERROR: could not find name key')
                coldata.append([True,str(b[use_name]).strip(),b['ra'],b['dec']])
            else:
                coldata.append([False, '--', np.nan, np.nan])

        coldata = list(map(list, zip(*coldata)))
        for i,col in enumerate(['_mask','','_ra','_dec']):
            table.add_column(Column(coldata[i], name=procname+col))

    table = sanitize_table(table)

    if kwargs['candidate']:
        table.write(kwargs['candidate'], format='ascii', overwrite=True)

    return(table)

def check_time(table, **kwargs):

    if kwargs['alert_time']:
        if 'time_mask' in table.keys():
            table.remove_column('time_mask')
        mask = []
        for i,row in enumerate(copy.copy(table)):
            t = Time(row['discovery_date'])

            rel_time = t.mjd - kwargs['alert_time'].mjd
            trange = kwargs['time']
            tlow=trange[0].to_value('day') ; thigh=trange[1].to_value('day')

            if (rel_time < tlow or rel_time > thigh):
                mask.append(True)
            else:
                mask.append(False)

    table.add_column(Column(mask, name='time_mask'))

    return(table)

def check_prob(table, **kwargs):

    if '2d_probability' in table.keys():
        if 'prob_mask' in table.keys():
            table.remove_column('prob_mask')
        prob = kwargs['probability']
        table.add_column(Column(table['2d_probability']>prob), name='prob_mask')

    return(table)

def check_photometry(table, **kwargs):

    # Import external and internal photometry
    photdata = None
    if 'photyse' in kwargs.keys() and kwargs['photyse']:
        #print('DOWNLOADING YSE \n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
        photdata = download_yse(kwargs['photyse'])
    if ('internal' in kwargs.keys() and kwargs['internal'] and
        os.path.exists(kwargs['internal'])):
        internal_photometry = ascii.read(kwargs['internal'])

        # Reformatting
        dates = [' '.join(Time(r['mjd'], format='mjd').fits.split('T'))
            for r in internal_photometry]
        internal_photometry.remove_column('mjd')
        internal_photometry.add_column(Column(dates, name='discovery_date'))
        internal_photometry.add_column(Column(dates, name='obsdate'))
        coords = SkyCoord(internal_photometry['ra'], internal_photometry['dec'],
            unit=(u.hour, u.deg))
        racol = Column([c.ra.degree for c in coords], name='ra')
        deccol = Column([c.dec.degree for c in coords], name='dec')
        internal_photometry['ra'] = racol
        internal_photometry['dec'] = deccol

        # Stack internal and external photometry
        if photdata:
            photdata = vstack([photdata, internal_photometry])
        else:
            photdata = internal_photometry

    # Reformat filters
    if photdata and 'filter' in photdata.keys():
        filts = []
        for row in photdata:
            filts.append(row['filter'].replace('-ZTF',''))
        photdata['filter'] = Column(filts, name='filter')

        for i,row in enumerate(photdata):
            if row['filter'].strip()=='Unknown':
                photdata['filter'][i]='r'

    if photdata:
        table = add_phot_data(table, photdata, **kwargs)
        table = add_phot_masks(table, 'var', **kwargs)
        table = add_phot_masks(table, 'absmag', **kwargs)

        mask = []
        for row in table:
            var = False
            for filt in table.meta['unique_filters']:
                var = var or row[filt+'_absmag_mask'] or row[filt+'_var_mask']

            mask.append(var)

        if 'photometry_mask' in table.keys():
            table.remove_column('photometry_mask')

        table.add_column(Column(mask, name='photometry_mask'))

    return(table)

def add_phot_masks(table, typ, **kwargs):
    var_filt = ''
    delta_mag = ''
    filts = table.meta['unique_filters']
    for filt in filts:
        mask = []
        for row in table:
            flag = False
            if not np.isnan(row[filt+'_'+typ]):
                if typ=='var':
                    flag=float(row[filt+'_var'])+float(row[filt+'_var_err'])<0.1
                if typ=='absmag':
                    flag=float(row[filt+'_absmag'])-\
                        0.1*float(row[filt+'_absmag_date'])+\
                        float(row[filt+'_absmag_err']) < -18.0

            mask.append(flag)

        if filt+'_'+typ+'_mask' in table.keys():
            table.remove_column(filt+'_var_mask')

        table.add_column(Column(mask, name=filt+'_'+typ+'_mask'))

    mask = [] ; use_filts = []
    for row in table:
        masks = [row[filt+'_'+typ+'_mask'] for filt in filts]
        if any(masks):
            mask.append(True)
            use_filts.append(filts[np.where(masks)[0][0]])
        else:
            mask.append(False)
            use_filts.append('--')

    if typ+'_mask' in table.keys():
        table.remove_column(typ+'_mask')
    if typ+'_filt' in table.keys():
        table.remove_column(typ+'_filt')

    table.add_column(Column(mask, name=typ+'_mask'))
    table.add_column(Column(use_filts, name=typ+'_filt'))

    return(table)

def step_message(step_num, proc, kwargs):
    msg_fmt = 'STEP {0: <3}: {1: <20} {2}'

    if step_num==0:
        msg = 'there are {curr:<4} viable candidates'
    else:
        msg =  'there are {curr:<4} viable candidates ({rem:<4} removed total;'
        msg += ' {this: <4} this step; {only: <4} only this step)'

    procname = format_proc(proc, upper=True)

    args = [str(step_num).zfill(3), procname, msg]
    m = msg_fmt.format(*args)
    print(m.format(**kwargs))

def get_kwargs(table, masks=[]):

    kwargs = {'curr': len(table)}

    mcopy = copy.copy(masks)

    if len(masks)>0:
        mcopy = np.array([np.array([bool(f) for f in m]) for m in mcopy])
        mask = np.any(np.array(mcopy), axis=0)
        this_mask = np.array(mcopy)[-1,:]
        if len(masks)>1:
            pre_mask = np.any(np.array(mcopy)[:-1,:], axis=0)
        else:
            pre_mask = np.array([False]*len(table))
        kwargs['curr']=len(table[~mask])
        kwargs['rem']=len(table[mask])
        kwargs['this']=len(table[this_mask & ~pre_mask])
        kwargs['only']=len(table[this_mask])

    return(kwargs)

def get_distance_modulus(row, avoid_mpc=True, typ='spec'):

    if (avoid_mpc and row['mpc_mask']) or typ=='--':
        if 'dl' in row.colnames and 'dl_err' in row.colnames:
            # dl and dl_err in units of Mpc
            dm = 5*np.log10(row['dl'])+25.0
            dm_err = 2.171472 * row['dl_err']/row['dl']
            return(dm, dm_err)

    if '_z' not in typ:
        typ=typ+'_z'

    if typ in row.colnames and typ+'_err' in row.colnames:
        try:
            z = float(row[typ])
            z_err = float(row[typ+'_err'])
            distance = cosmo.luminosity_distance(z).to_value('parsec')
            distance_err = distance * z_err / z

            dm = 5*np.log10(distance)-5.0
            dm_err = 2.171472 * distance_err / distance

            dm_err = np.sqrt(0.1**2+dm_err**2)

            return(dm, dm_err)

        except ValueError:
            return(np.nan, np.nan)

    if 'dl' in row.colnames and 'dl_err' in row.colnames:
        dm = 5*np.log10(row['dl'])+25.0
        dm_err = 2.171472 * row['dl_err']/row['dl']
        return(dm, dm_err)
    else:
        return(np.nan, np.nan)

def pick_best_redshift(table, **kwargs):

    redshift_cols = []
    for key in table.keys():
        if key.endswith('_z') and key+'_err' in table.keys():
            redshift_cols.append(key.replace('_z',''))

    use_redshift = ['--']*len(table)

    for i,row in enumerate(table):

        for col in redshift_cols:
            try:
                z = float(row[col+'_z'])
                z_err = float(row[col+'_z_err'])

                if np.isnan(z) or np.isnan(z_err):
                    continue

                coord = parse_coord(row['ra'], row['dec'])
            except ValueError:
                continue

            if use_redshift[i]=='spec' or col=='spec':
                use_redshift[i]='spec'
                continue

            if use_redshift[i]=='--':
                use_redshift[i] = col
            else:
                comp_z = float(row[use_redshift[i]+'_z'])
                comp_z_err = float(row[use_redshift[i]+'_z_err'])

                comp_ra = float(row[use_redshift[i]+'_ra'])
                comp_dec = float(row[use_redshift[i]+'_dec'])
                comp_coord = parse_coord(comp_ra, comp_dec)
                comp_sep = comp_coord.separation(coord).to_value('radian')
                comp_proj_sep = comp_sep * cosmo.luminosity_distance(comp_z)

                z_ra = float(row[col+'_ra'])
                z_dec = float(row[col+'_ra'])
                z_coord = parse_coord(z_ra, z_dec)

                if z_coord is None:
                    continue

                z_sep = z_coord.separation(coord).to_value('radian')
                z_proj_sep = z_sep * cosmo.luminosity_distance(z)

                # What if the sources are the same?  Use one with smaller error
                if z_coord.separation(comp_coord).to_value('arcsec') < 2.0:
                    if z_err / z < comp_z_err / comp_z:
                        use_redshift[i] = col
                # Otherwise use the one with the smaller projected separation
                else:
                    if z_proj_sep < comp_proj_sep:
                        use_redshift[i] = col

    if kwargs['redshift']['use_name'] in table.keys():
        table.remove_column(kwargs['redshift']['use_name'])
    if 'redshift_mask' in table.keys():
        table.remove_column('redshift_mask')

    # Populate the mask criterion from the other redshift masks
    cut_redshift = []
    for i,row in enumerate(table):
        if use_redshift[i]!='--' and row[use_redshift[i]+'_mask']:
            cut_redshift.append(True)
        else:
            cut_redshift.append(False)

    table.add_column(Column(use_redshift, name=kwargs['redshift']['use_name']))
    table.add_column(Column(cut_redshift, name='redshift_mask'))

    return(table)

def get_rf(filt):
    rf = extCurves.avebvsf
    filt = filt.strip()
    if filt in ['u','g','r','i','z','up','gp','rp','ip','ip','zp']:
        filt = filt.replace('p','')
        return(rf['SDSS '+filt])
    if filt=='y':
        return(rf['PS1 y'])
    if filt in ['U','B','V','R','I']:
        return(rf['CTIO '+filt])
    elif filt=='w':
        return(0.051/0.058)
    else:
        try:
            red = rf[filt]
            return(red)
        except:
            return(None)

def add_phot_data(table, photdata, output_phot='photlist.dat', **kwargs):

    if 'dm' in table.keys():
        table.remove_column('dm')
    if 'dm_err' in table.keys():
        table.remove_column('dm_err')

    distance_modulus = [] ; distance_modulus_err = []
    for row in table:
        dm, dmerr = get_distance_modulus(row, typ=row['use_redshift'].lower())
        distance_modulus.append(dm) ; distance_modulus_err.append(dmerr)

    # Substitute rp->r, gp->g, ip->i, up->u
    mask = photdata['filter']=='rp'
    photdata[mask]['filter']='r'
    mask = photdata['filter']=='ip'
    photdata[mask]['filter']='i'
    mask = photdata['filter']=='gp'
    photdata[mask]['filter']='g'
    mask = photdata['filter']=='up'
    photdata[mask]['filter']='u'

    table.add_column(Column(distance_modulus, name='dm'))
    table.add_column(Column(distance_modulus_err, name='dm_err'))

    if output_phot:
        output_phot_file = open(output_phot, 'w')

    if ('mw_ebv' not in table.keys() or
        any([np.isnan(val) for val in table['mw_ebv']])):
        mw_ebv = []
        for row in table:
            try:
                val = sfd(parse_coord(row['ra'], row['dec']))
                mw_ebv.append(val)
            except:
                mw_ebv.append(np.nan)
        table.add_column(Column(mw_ebv, name='mw_ebv'))

    for filt in kwargs['avoid_filters']:
        mask = photdata['filter']!=filt
        photdata = photdata[mask]

    unique_filters = np.unique(photdata['filter'].data)
    table.meta['unique_filters']=np.array(unique_filters)

    NFILTS = len(unique_filters)
    NTRANSIENTS = len(table)
    size = (NTRANSIENTS, NFILTS)

    absmag = np.empty(size)
    absmag_date = np.empty(size)
    absmag_err = np.empty(size)
    variability = np.empty(size)
    variability_err = np.empty(size)

    for i,row in enumerate(table):
        mask = ((photdata['name']==row['name']) &
                (photdata['magerr']<1.0))
        photometry = photdata[mask]

        band_id_dict={'119': 'DECam','120': 'DECam','40': 'PS1',
                        '39': 'PS1','5': 'Swope','109': 'Thacher',
                        '116': 'LCO','10': 'Nickel','9': 'Nickel',
                        '115': 'LCO', '113': 'LCO', '78': 'ZTF',
                        '6': 'Swope', '3': 'Swope', '4': 'Swope',
                        '--': 'Other','77':'ZTF'}

        for data_point in photometry:
            if np.isnan(float(data_point['magerr'])):
                        continue
            if (not np.isnan(data_point['magerr'])
                and 'nan' not in str(data_point['magerr'])
                and Time(data_point['obsdate']).mjd>58709.00
                and str(data_point['band_id'])!='11'
                and str(data_point['band_id'])!='118'
                and str(data_point['band_id'])!='117'
                and Time(data_point['obsdate']).mjd<58757.00):
                pfilt = data_point['filter'].replace('p','')
                band_name = '--'
                if str(data_point['band_id']) in band_id_dict.keys():
                    band_name=band_id_dict[str(data_point['band_id']).strip()]
                name = row['name']
                mjd_str = '%5.4f'%Time(data_point['obsdate']).mjd
                mag_str = '%2.2f'%data_point['mag']
                magerr_str = '%1.2f'%data_point['magerr']
                band_id_str = str(data_point['band_id']).strip()

                outline = '{0} {1} {2} {3} {4} {5} {6}'
                outline = outline.format(name, mjd_str, pfilt,
                    mag_str, magerr_str, band_name, band_id_str)
                print(outline)
                if output_phot:
                    output_phot_file.write(outline+'\n')

        dm = row['dm'] ; dm_err = row['dm_err'] ; ebv = row['mw_ebv']

        for j,filt in enumerate(unique_filters):
            mask = np.array([filt.strip() in photometry['filter']])

            if kwargs['alert_time']:
                a = kwargs['alert_time']
                times = np.array([Time(r['obsdate']).mjd-Time(a).mjd
                        for r in photometry])
                mask = mask & (times > 0)

            filter_photometry = photometry[mask]

            if len(filter_photometry)>0:

                filter_photometry.sort('obsdate')
                initial_data = filter_photometry[0]

                mag = initial_data['mag']
                magerr = initial_data['magerr']

                if np.isnan(magerr) or str(magerr)=='--':
                    magerr = 0.2

                if kwargs['alert_time']:
                    obsdate = Time(initial_data['obsdate']).mjd - a.mjd
                else:
                    obsdate = Time(initial_data['obsdate']).mjd

                rf = get_rf(filt)

                absmag[i,j]=mag-dm-rf*ebv
                absmag_date[i,j]=obsdate
                absmag_err[i,j]=np.sqrt(magerr**2 + dmerr**2)

                delta_mag = np.nan ; unc = np.nan
                obsdates = filter_photometry['obsdate']
                times = np.array([Time(r).mjd for r in obsdates])
                if kwargs['alert_time']: times = times - a.mjd

                if np.max(times)-np.min(times) > 0.1:
                    mags = np.array(filter_photometry['mag'].data)
                    errs = np.array(filter_photometry['magerr'].data)
                    popt, pcov = curve_fit(lin_decline, times, mags, sigma=errs)

                    unc = pcov[1,1]+0.05
                    if np.isinf(unc):
                        unc = np.sqrt(np.sum([e**2 for e in errs]))

                    delta_mag = popt[1]

                variability[i,j]=delta_mag
                variability_err[i,j]=unc

            else:
                absmag[i,j]=np.nan
                absmag_date[i,j]=np.nan
                absmag_err[i,j]=np.nan
                variability[i,j]=np.nan
                variability_err[i,j]=np.nan

    for j,filt in enumerate(unique_filters):

        for key in ['_absmag','_absmag_err','_absmag_date','_var','_var_err']:
            if filt+key in table.keys():
                table.remove_column(filt+key)

        table.add_column(Column(absmag[:,j], name=filt+'_absmag'))
        table.add_column(Column(absmag_err[:,j], name=filt+'_absmag_err'))
        table.add_column(Column(absmag_date[:,j], name=filt+'_absmag_date'))
        table.add_column(Column(variability[:,j], name=filt+'_var'))
        table.add_column(Column(variability_err[:,j], name=filt+'_var_err'))

    table = sanitize_table(table)

    if output_phot:
        output_phot_file.close()

    if kwargs['candidate']:
        table.write(kwargs['candidate'], format='ascii', overwrite=True)

    return(table)
