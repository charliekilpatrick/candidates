parse_candidates.py transients_list_220612_220613.txt --no-lc-error
# On ziggy in scratch directory
python ../analysis.py candidates.csv --redo --steps import ps1dr2 astcheck gaia --candidate-format
