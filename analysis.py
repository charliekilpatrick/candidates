import warnings, os, sys, copy, numpy as np
warnings.filterwarnings('ignore')

import utilities.utilities as util
import options

if __name__=="__main__":

    args = options.add_options()
    kwargs = options.parse_kwargs(args)
    steps = options.parse_steps(args)

    table = None ; masks = []
    finished_steps = []
    for i,proc in enumerate(steps):
        table = util.add_data(table, proc, **kwargs)

        if i>0: masks.append(table[util.format_proc(proc)+'_mask'].data)
        if len(masks)>0:
            mcopy = copy.copy(masks)
            mcopy = np.array([np.array([bool(f) for f in m]) for m in mcopy])
            mask = np.any(np.array(mcopy), axis=0)
            subtable=table[~mask]

            outdatafile = os.path.join(kwargs['event_dir'],
                util.format_proc(proc)+'_table.dat')
            subtable.write(outdatafile, format='ascii', overwrite=True)

        util.step_message(i, proc, util.get_kwargs(table, masks=masks))
        finished_steps.append(util.format_proc(proc))

        if 'time' in finished_steps and 'prob' in finished_steps:
            if 'initmask' in table.meta.keys() and table.meta['initmask']:
                continue
            else:
                tmask = (~table['time_mask']) & (~table['prob_mask'])
                table = table[tmask]
                table.meta['initmask']=True

    table.meta['use_masks']=masks
    if args.latex:
        util.output_latex_table(table, **kwargs)
    if args.candidate_format:
        util.output_candidate_table(table, **kwargs)
