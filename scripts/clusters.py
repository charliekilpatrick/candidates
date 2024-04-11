#!/usr/bin/env python

import numpy as np
import os
import sys
import warnings
import glob

from astropy.table import Table, vstack, Column
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt

photcode_dict={'01':'U',
               '02':'B',
               '03':'V',
               '04':'R',
               '05':'I',
               '12':'u',
               '13':'g',
               '14':'r',
               '15':'i',
               '16':'z'}

def add_options():
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("objname", type=str,
        help="Name of source at the start of objname to glob cluster files.")
    parser.add_argument("ra", type=str,
        help="RA for source to search for in clusters files.")
    parser.add_argument("dec", type=str,
        help="Dec for source to search for in clusters files.")
    parser.add_argument("--plot", default=False, action='store_true',
        help="Output a figure for the parsed photometry.")
    parser.add_argument("--input", default=None, type=str,
        help="File name of an input file to use instead of parsing data to plot.")
    parser.add_argument("--mjd-range", nargs=2,
        help="Restrict plot MJD range to these values.")
    parser.add_argument("--filters", nargs='+',
        help="Restrict plot to these filters.")
    parser.add_argument("--compare", default=None, type=str,
        help="A file path for a photometry file to compare to given data in plot.")
    parser.add_argument("--mag-range", nargs=2,
        help="Restrict plot magnitude range to these values.")

    args = parser.parse_args()

    return(args)

def parse_coord(ra, dec):
    if ':' in ra and ':' in dec:
        unit=(u.hour, u.deg)
    else:
        unit=(u.deg, u.deg)

    try:
        coord = SkyCoord(ra, dec, unit=unit, frame='icrs')
        return(coord)
    except:
        return(None)

def transpose(data):
    return(list(map(list, zip(*data))))

def import_cluster_file(filename):
        idx=0
        table_count=0
        filedata=[]
        photdata=[]
        canddata=[]
        phottabs={}
        with open(filename,'r') as f:
            for line in f:
                if line.startswith('#'):
                    if table_count==0:
                        filehead=line.replace('#','').split()
                    elif table_count==1:
                        candhead=line.replace('#','').split()
                    else:
                        phothead=line.replace('#','').split()
                    
                    if photdata:
                        photkey = str(photdata[0][0].strip())
                        if photkey in phottabs.keys():
                            phottabs[photkey].extend(photdata)
                        else:
                            phottabs[photkey]=photdata
                        photdata=[]
                    
                    table_count += 1
                    continue
                if table_count==1:
                    filedata.append(line.split())
                elif table_count==2:
                    canddata.append(line.split())
                else:
                    photdata.append(line.split())
            
            if photdata:
                photkey = str(photdata[0][0].strip())
                if photkey in phottabs.keys():
                    phottabs[photkey].extend(photdata)
                else:
                    phottabs[photkey]=photdata
                photdata=[]

        filetable = Table(transpose(filedata), names=filehead)
        candtable = Table(transpose(canddata), names=candhead)
        for key in phottabs.keys():
            table = Table(transpose(phottabs[key]), names=phothead)
            for subkey in table.keys(): table.rename_column(subkey, subkey.lower())
            phottabs[key] = table

        return(filetable, candtable, phottabs)

def get_all_tables(table_dir):
    files = None
    candidates = None

    globstr=table_dir+'*.clusters'
    if len(table_dir+'*.clusters')==0:
        raise Exception(f'ERROR: no clusters files following {globstr}')

    for file in glob.glob(globstr):
        filetable, candtable, phottabs = import_cluster_file(file)
        
        if not files:
            files = Table(filetable)
        else:
            files = vstack([files, filetable])
        
        if not candidates:
            candidates = Table(candtable)
        else:
            candidates = vstack([candidates, candtable])

    for key in files.keys(): files.rename_column(key, key.lower())
    for key in candidates.keys(): candidates.rename_column(key, key.lower())

    return(files, candidates, phottabs)

def identify_candidate(ctable, ra, dec):
    coord = parse_coord(ra, dec)

    coords = np.array([parse_coord(row['raaverage'], row['decaverage'])
        for row in ctable])

    seps = np.array([coord.separation(c).to_value(u.arcsec) for c in coords])
    mask = seps < 2.0

    if len(ctable[mask])>0:
        subtable = ctable[mask]

        coords = np.array([parse_coord(row['raaverage'], row['decaverage'])
            for row in subtable])
        seps = np.array([coord.separation(c).to_value(u.arcsec) for c in coords])

        sep = Column(seps, name='separation')
        subtable.add_column(sep)
        subtable.sort('separation')

        return(subtable[0]['id'])
    else:
        return(None)

def get_photometry(ftable, ptable, idnum, forced=False):

    subtable = ptable[idnum]

    if forced: tmask = subtable['type']=='0x00000011'
    else: tmask = subtable['type']=='0x00000001'

    subtable = subtable[tmask]

    newtable = Table([[0.],['X'*100],[0.],[0.]],
        names=('mjd','filter','mag','magerr')).copy()[:0]

    for row in subtable:

        mag = float(row['m'])
        merr = float(row['dm'])
        ferr = float(row['dflux'])

        fmask = ftable['cmpfile']==row['cmpfile']
        if ftable[fmask][0]['convol00'].lower() not in ['image','template']:
            continue

        zpt = float(ftable[fmask][0]['zptmagav'])
        zerr = float(ftable[fmask][0]['tzptmuce'])

        mjd = ftable[fmask][0]['mjd']
        photcode = ftable[fmask][0]['photcode']

        filt = photcode_dict[photcode[-2:]]

        mag = mag + zpt
        merr = np.sqrt(merr**2 + zerr**2)

        if merr>0.33:
            mag = -2.5 * np.log10(3 * ferr) + zpt
            merr = 0.0

        mag = '%2.4f'%float(mag)
        merr = '%2.4f'%float(merr)

        print(mjd,filt,mag,merr)

        newtable.add_row([mjd,filt,mag,merr])

    return(newtable)

def plot_data(outdata, objname, set_mjd_range=None, set_mag_range=None,
    set_filters=None, compare=None):

    colors = {'u': 'violet','g': 'green','r': 'red','i':'orange', 'z':'black',
        'cyan-ATLAS':'cyan','g-ZTF':'cyan','r-ZTF':'red','orange-ATLAS':'brown',
        'B': 'blue', 'V': 'darkgreen', 'R': 'darkred','U':'magenta',
        'zs':'black'}

    good_mag = (outdata['magerr']!=0.0) & (outdata['magerr']<0.33)

    mjd_range = np.max(outdata['mjd'])-np.min(outdata['mjd'])
    mag_range = np.max(outdata['mag'][good_mag])-np.min(outdata['mag'][good_mag])
    mjd_buffer = 0.05 * mjd_range
    mag_buffer = 0.05 * mag_range

    fig, ax = plt.subplots()
    for filt in np.unique(outdata['filter']):

        if set_filters is not None:
            if filt not in set_filters:
                continue

        mask = outdata['filter']==filt
        label = filt
        for row in outdata[mask]:
            if row['magerr']>0:
                marker='o'
            else:
                marker='v'
            ax.errorbar(row['mjd'], row['mag'], yerr=row['magerr'], marker=marker,
                color=colors[filt], label=label, markeredgecolor='k')
            if label:
                label = None

    if compare:
        compare = Table.read(compare, format='ascii', names=('mjd','filter',
            'mag','magerr'))
        for filt in np.unique(compare['filter']):
            
            if set_filters is not None:
                if filt not in set_filters:
                    continue

            mask = compare['filter']==filt
            label = filt
            for row in compare[mask]:
                if row['magerr']>0:
                    marker='s'
                else:
                    marker='v'
                ax.errorbar(row['mjd'], row['mag'], yerr=row['magerr'], marker=marker,
                    color=colors[filt], label=label, markeredgecolor='k')
                if label:
                    label = None

    ax.set_ylim([np.max(outdata['mag'])+mag_buffer, np.min(outdata['mag'])-mag_buffer])
    ax.set_xlim([np.min(outdata['mjd'])-mjd_buffer, np.max(outdata['mjd'])+mjd_buffer])

    if set_mjd_range is not None:
        set_mjd_range = [float(f) for f in set_mjd_range]
        ax.set_xlim(set_mjd_range)

    if set_mag_range is not None:
        set_mag_range = [float(f) for f in set_mag_range]
        ax.set_ylim(set_mag_range)

    plt.legend()

    plt.tight_layout()
    plt.savefig(os.path.basename(objname)+'.pdf')

def reformat_outdata(outdata):

    for i,row in enumerate(outdata):

        filt = row['filter']
        if filt.endswith('p'):
            newfilt = filt.replace('p','')
            outdata[i]['filter']=newfilt

    return(outdata)

if __name__=="__main__":

    if len(sys.argv) < 4: args = add_options()

    coord = parse_coord(sys.argv[2], sys.argv[3])
    sys.argv[2] = str(coord.ra.degree) ; sys.argv[3] = str(coord.dec.degree)
    args = add_options()

    if args.input is None:

        files, candidates, photometry = get_all_tables(args.objname)

        idnum=identify_candidate(candidates, args.ra, args.dec)

        outdata = get_photometry(files, photometry, idnum, forced=True)
    else:
        outdata = Table.read(args.input, format='ascii',
            names=('mjd','filter','mag','magerr'))

    outdata = reformat_outdata(outdata)

    if args.plot:
        plot_data(outdata, args.objname, set_mjd_range=args.mjd_range,
            set_mag_range=args.mag_range,
            set_filters=args.filters, compare=args.compare)

    with open(os.path.basename(args.objname)+'.phot','w') as f:
        for filt in sorted(np.unique(outdata['filter'])):
            subdata = outdata[outdata['filter']==filt]
            subdata.sort('mjd')
            for row in subdata:
                output = '%5.4f %s %2.4f %.4f \n'%(row['mjd'],row['filter'],row['mag'],row['magerr'])
                f.write(output)




