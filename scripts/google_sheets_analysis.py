import pickle
import os
import sys
import requests
import time
import string

from requests.auth import HTTPBasicAuth

from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
from dateutil.parser import parse

from astropy.io import ascii
from astropy.table import Table, vstack, unique
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u

from collections import Counter
import numpy as np

import hashlib

from gw190425_sheets import *

def transpose(l):
    return(list(map(list, zip(*l))))

def initiate_gsheet(token_file):
    creds = None
    if os.path.exists(token_file):
        with open(token_file, 'rb') as token:
                creds = pickle.load(token)

    service = build('sheets', 'v4', credentials=creds)
    sheet = service.spreadsheets()

    return(sheet)

def is_number(val):
    try:
        x = float(val)
        return(True)
    except ValueError:
        return(False)


def get_websniff_data(sheet, spreadsheetId):

    for i in np.arange(3):
        try:
            response = sheet.values().get(spreadsheetId=spreadsheetId, 
                range='WEBSNIFF').execute()
            break
        except googleapiclient.errors.HttpError:
            print(f'WARNING: API limit exceeded.  Waiting 1 minute...')
            time.sleep(60)

    table = Table([['X'*100],[0],['X'*100],['X'*100],['X'*100]],
        names=('page','ncandidates','sniff1','sniff2','sniff3')).copy()[:0]

    if 'values' not in response.keys():
        return(table)

    for data in response['values']:

        if len(data)<3: continue

        fieldname = data[0]
        num = fieldname.split('.')[-1]

        if '.' in fieldname and '_' in fieldname and is_number(num):
            row_data = [fieldname]
            row_data.append(int(data[2]))
            if len(data)>3:
                row_data.append(data[3])
            if len(data)>4:
                row_data.append(data[4])
            if len(data)>5:
                row_data.append(data[5])

            while len(row_data)<5: row_data.append('')

            table.add_row(row_data)

    return(table)

def get_candidates(sheet, spreadsheetId, name):

    for i in np.arange(3):
        try:
            response = sheet.values().get(spreadsheetId=spreadsheetId, 
                range=name).execute()
            break
        except googleapiclient.errors.HttpError:
            print(f'WARNING: API limit exceeded.  Waiting 1 minute...')
            time.sleep(60)

    candidates=[]

    if 'values' not in response.keys():
        return([])

    for data in response['values']:
        if len(data)==0: continue
        if 'http://' in data[0]:
            candidates.append(data[0])

    return(candidates)

def parse_names(table):

    unique_names=[]
    for row in table:
        for key in ['sniff1','sniff2','sniff3']:
            name = row[key].split()[0].strip()
            if not name: continue
            if name not in unique_names:
                unique_names.append(name)

    return(unique_names)

def parse_telescope_candidates(sheet, spreadsheetId):

    table = get_websniff_data(sheet, spreadsheetId)
    names = parse_names(table)

    candidates={}
    for name in names:
        candidates[name]=get_candidates(sheet, spreadsheetId, name)

    good_candidates={}
    all_candidates = []
    for row in table:
        page=row['page']
        if not row['sniff1'] or not row['sniff2'] or not row['sniff3']:
            print(f'{page} has less than 3 people who sniffed!')
        else:
            page_candidates = []
            good_candidates[page]=Table([['X'*100],[0]],
                names=('candidate','count')).copy()[:0]
            for sniffer in ['sniff1','sniff2','sniff3']:
                # Get unique so we don't have duplicates from a single sniffer
                name = row[sniffer].split()[0].strip()
                sniff_candidates=np.unique(
                    [c for c in candidates[name] if page in c])
                page_candidates.extend(sniff_candidates)

            # Now we can simply count the incidence of each candidate across, 
            # since duplicate entries will correspond to different people 
            # sniffing it
            c = Counter(page_candidates)
            for key,val in zip(c.keys(), c.values()):
                good_candidates[page].add_row([key, val])

            good_candidates[page].sort('count')

            # Initial masking out of candidates - gets rid of duplicates if we 
            # accidentally added them more than once
            if len(good_candidates[page])>0:
                good_candidates[page] = unique(good_candidates[page],
                    keys=['candidate'], keep='last')

            mask = good_candidates[page]['count']>1
            for row in good_candidates[page][mask]:
                # Masks out multiple candidates a second time, e.g., for google
                # sheets pages where websniff pages are repeated for
                # large (>100) numbers of candidates
                if row['candidate'] not in [a[1] for a in all_candidates]:
                    all_candidates.append([page, row['candidate']])

    if len(all_candidates)==0:
        table = Table([['X'*100],['X'*100]],
            names=('page', 'candidate')).copy()[:0]
    else:
        table = Table(transpose(all_candidates), names=('page', 'candidate'))
    
    return(table)

def get_lightcurve_file(candidate_url, typ='unforced', verbose=True):

    base_url = candidate_url.split('.html')[0]
    cand_id = candidate_url.split('#')[-1]

    lc_url = f'{base_url}_cand{cand_id}.{typ}.difflc.txt'
    lc_url = lc_url.replace(' ','')

    if 'http://' in lc_url: lc_url=lc_url.replace('http://','https://')

    if verbose: print(f'Trying to get: {lc_url}')

    if 'ziggy' in lc_url:
        auth = HTTPBasicAuth(os.environ['ZIGGY_USER'], 
            os.environ['ZIGGY_PASSWORD'])
    else:
        auth = None

    r = requests.get(lc_url, auth=auth)

    if verbose: print('Status code:',r.status_code)

    if r.status_code==200:
        table = ascii.read(r.text)
        return(table)
    else:
        return(None)


def parse_name(candidate, mjd, ra, dec):

    candidate = candidate.lower()
    cand_name = ''
    if 'nickel' in candidate:
        cand_name = 'NGW'
    elif 'swope' in candidate:
        cand_name = 'SSS'
    elif 'thacher' in candidate:
        cand_name = 'TGW'
    elif 'andicam' in candidate:
        cand_name = 'AGW'

    date_str = Time(mjd, format='mjd').datetime.strftime('%Y')[2:]
    cand_name = cand_name+date_str

    mystring=ra+dec
    res = bytes(mystring, 'utf-8')
    rand = str(int(hashlib.sha256(res).hexdigest(), base=16) % 2345679823456798)

    letters = string.ascii_lowercase
    rand_letters = ''
    for i in np.arange(8):
        val = rand[2*i:2*i+2]
        try:
            val = int(val) % 26
        except ValueError:
            val = 0
        rand_letters += letters[val]

    cand_name = cand_name + rand_letters

    return(cand_name)

sheet = initiate_gsheet(os.environ['GSHEET_TOKEN_FILE'])

table = None
for key in all_sheets.keys():
    new = parse_telescope_candidates(sheet, all_sheets[key])
    if not table:
        table = new
    else:
        table = vstack([table, new])

out_candidate_table = Table([['X'*100],[0.],[0.],[0.],['X'*10],[0.],[0.],
    ['X'*100]], names=('name','mjd','ra','dec','filter','mag','mag_err',
    'candidate_url')).copy()[:0]
for row in table:
    lc = get_lightcurve_file(row['candidate'], verbose=True)

    if lc is None:
        cand=row['candidate']
        raise Exception(f'ERROR: could not download lc file for {cand}')

    lc.sort('MJD')
    lc = lc[0]
    name = parse_name(row['candidate'], lc['MJD'], lc['ra'], lc['dec'])

    coord = SkyCoord(lc['ra'], lc['dec'], unit=(u.hour, u.deg))

    if '-' in str(lc['m']):
        mag = 0.0
    else:
        mag = float(lc['m'])

    if '-' in str(lc['dm']):
        dmag = 0.0
    else:
        dmag = float(lc['dm'])

    out_candidate_table.add_row([name, lc['MJD'], coord.ra.degree,
        coord.dec.degree, lc['filt'], mag, dmag, row['candidate']])

out_candidate_table.write('gw190425_candidates.csv', format='csv',
    overwrite=True)
