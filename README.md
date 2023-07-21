# candidates

Methods for vetting candidate optical transients as extragalactic and associated with gravitational wave (GW) events.  Includes tools for comparing to bayestar GW maps within the localization region, occurrence from detection of the event, coincidence with a minor planet, Gaia or Pan-STARRS point source, coincidence with ASAS-SN variability, association with a host with a NED, PS1-STRM, Legacy, or 2MASS redshift, and photometric evolution.

## Installation

Provided is a `requirements.txt` and `candidates.yml` file for installation of dependencies via pip or conda.  These have been tested on conda v23.5.2 and Ubuntu v20.04.4.

## Options

```
usage: analysis.py [-h] [--gw-event GW_EVENT] [--cut-time CUT_TIME] [-r] [-s SHIBBOLETH] [-m MAP_TYPE] [--steps STEPS [STEPS ...]]
                   [--max-date MAX_DATE] [--max-prob MAX_PROB] [--distance-distribution DISTANCE_DISTRIBUTION] [-z REDSHIFT REDSHIFT]
                   [--mpc-radius MPC_RADIUS] [--astcheck-radius ASTCHECK_RADIUS] [--gaia-radius GAIA_RADIUS] [--asassn-radius ASASSN_RADIUS]
                   [--ps1-radius PS1_RADIUS] [--tns-radius TNS_RADIUS] [--yse-radius YSE_RADIUS] [--gal-radius GAL_RADIUS]
                   [--max-gal-radius MAX_GAL_RADIUS] [--phot-bright PHOT_BRIGHT] [--phot-min-decline PHOT_MIN_DECLINE]
                   [--phot-color PHOT_COLOR] [--yse-phot YSE_PHOT] [--redshift-methods REDSHIFT_METHODS [REDSHIFT_METHODS ...]]
                   [--tns-api-key TNS_API_KEY] [--tns-bot-name TNS_BOT_NAME] [--tns-bot-id TNS_BOT_ID] [--latex] [--candidate-format]
                   [--obscode OBSCODE] [--verbose]
                   sources

positional arguments:
  sources               Source of transient/candidate data to analyze. Can be a path to file or a URL.

options:
  -h, --help            show this help message and exit
  --gw-event GW_EVENT   Input event name. Must be parseable from GraceDB (required for probability cut).
  --cut-time CUT_TIME   Override event time with input value (usually parsed from Grace DB).
  -r, --redo            Redo all steps if they have not yet been run.
  -s SHIBBOLETH, --shibboleth SHIBBOLETH
                        Shibboleth file for sub-routines that require authentication.
  -m MAP_TYPE, --map-type MAP_TYPE
                        Map type for analysis (bayestar should always exist).
  --steps STEPS [STEPS ...]
                        List of steps to check before validating a candidate.
  --max-date MAX_DATE   Maximum days from GW discovery to analyze candidates (in days).
  --max-prob MAX_PROB   Maximum cumulative probability (from 0-1) to consider for valid candidates.
  --distance-distribution DISTANCE_DISTRIBUTION
                        Distance distribution (in sigma; 0-infinity) to consider for candidate analysis.
  -z REDSHIFT REDSHIFT, --redshift REDSHIFT REDSHIFT
                        Redshift range in which to search for candidates.
  --mpc-radius MPC_RADIUS
                        Radius for minor planet checker (in arcsec).
  --astcheck-radius ASTCHECK_RADIUS
                        Radius for astcheck minor planet code (in arcsec).
  --gaia-radius GAIA_RADIUS
                        Radius for crossmatching to Gaia stars (in arcsec).
  --asassn-radius ASASSN_RADIUS
                        Radius for crossmatching to ASASSN variables (in arcsec).
  --ps1-radius PS1_RADIUS
                        Radius for crossmatching to PS1 DR2 stars (in arcsec).
  --tns-radius TNS_RADIUS
                        Radius for crossmatching to TNS transients (in arcsec).
  --yse-radius YSE_RADIUS
                        Radius for crossmatching to YSE transients (in arcsec).
  --gal-radius GAL_RADIUS
                        Radius for crossmatching to candidate host galaxies (in arcsec).
  --max-gal-radius MAX_GAL_RADIUS, --max-galaxy-radius MAX_GAL_RADIUS
                        Maximum separation from a galaxy with a known distance/redshift (to rule out chance coincidence; in parsec).
  --phot-bright PHOT_BRIGHT
                        For candidates with known distance, maximum apparent magnitude.
  --phot-min-decline PHOT_MIN_DECLINE
                        Minimum decline rate (mag/day) of candidates to be considered.
  --phot-color PHOT_COLOR
                        Minimum color (g-r mag) of candidates to be considered.
  --yse-phot YSE_PHOT   ID of YSE-PZ photometry SQL query. Default is not to use one.
  --redshift-methods REDSHIFT_METHODS [REDSHIFT_METHODS ...]
                        List of sources to check for host redshifts [spec,ps1strm,2mpz,legacy]
  --tns-api-key TNS_API_KEY
                        TNS API key - required for TNS queries.
  --tns-bot-name TNS_BOT_NAME
                        TNS API key - required for TNS queries.
  --tns-bot-id TNS_BOT_ID
                        TNS API key - required for TNS queries.
  --latex               Output latex-formatted candidate file.
  --candidate-format    Output candidate-formatted file.
  --obscode OBSCODE     Observatory code for astcheck functionality. For a list, see https://minorplanetcenter.net//iau/lists/ObsCodes.html.
  --verbose             Verbose output from candidate vetting methods.
```

## Contact

For all questions, comments, suggestions, and bugs related to this script, please contact Charlie Kilpatrick at ckilpatrick@northwestern.edu.
