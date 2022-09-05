import os

from datetime import datetime
import mdtraj
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

from pymd.mdanalysis.postprocess import PostProcess

def now():
    return datetime.now().strftime("%d/%m/%Y %H:%M:%S")

def logger(rep, job):
    log = '{}.{}.{}.log'.format(rep.parent.name, rep.id, job)
    log = os.path.join(rep.parent.scripts, log)
    f = open(log, 'w')
    f.write('Created on: {}\n'.format(now()))
    f.write('Created by: pymd.mdrun.{}\n\n'.format(job))
    return f


def pbc(rep):
    log = logger(rep, 'pbc')
    log.write('Loading trajectory... start time {}\n'.format(now()))
    try:
        traj = mdtraj.load(rep.xtc, top=rep.gro)
    except Exception as e:
        log.write('Loading trajectory failed.\n{}\n'.format(e))
        log.close()
        return e
    log.write('Trajectory loaded... end time {}\n\n'.format(now()))
    log.write('Imaging trajectory... start time {}\n'.format(now()))
    try:
        traj = traj.image_molecules(make_whole=True)
    except Exception as e:
        log.write('Imaging trajectory failed.\n{}\n'.format(e))
        log.close()
        return e
    log.write('Imaging complete... end time {}\n\n'.format(now()))
    if rep.parent.job_params['output'] is None:
        output = os.path.join(rep.root, 'pbc.xtc')
    else:
        output = os.patch.join(rep.root, rep.parent.job_params['output'])
    log.write('Writing corrected trajectory to {}... start time {}\n'.format(output, now()))
    try:
        traj.save_xtc(output)
    except Exception as e:
        log.write('Writing trajectory failed. \n{}\n'.format(e))
        log.close()
        return e
    log.write('Completed writing trajectory to {}.... end time {}\n'.format(output, now()))
    log.write('Process complete')
    log.close()
    return output

def rmsd(rep):
    # get job params 
    inp = os.path.join(rep.root, rep.parent.job_params['inp'])
    top = os.path.join(rep.root, rep.parent.job_params['top'])
    if top == 'system':
        top = os.path.join(rep.root, rep.parent.gro)
    precentered = rep.parent.job_params['precentered']
    selections = rep.parent.job_params['selections']
    output = os.path.join(rep.rmsd.root, rep.parent.job_params['output'])
    ref_index = rep.parent.job_params['reference_index']
    log = logger(rep, 'rmsd')
    log.write('Loading trajectory... start time {}\n'.format(now()))
    try:
        traj = mdtraj.load(inp, top=top)
    except Exception as e:
        log.write('Loading trajectory failed.\n{}\n'.format(e))
        log.close()
        return e
    log.write('Trajectory loaded... end time {}\n\n'.format(now()))
    log.write('RMSD\n')
    log.write('**********\n')
    log.write('Job Parameters:\n')
    log.write('Trajectory: {}\n'.format(inp))
    log.write('Topology: {}\n'.format(top))
    log.write('Selection(s): {}\n'.format(selections[0]))
    if len(selections) > 1:
        for sele in selections[1:]:
            log.write('\t\t{}\n'.format(sele))
    log.write('Reference Index: {}\n'.format(ref_index))
    log.write('Precentered Trajectory: {}\n'.format(precentered))
    log.write('Output: {}\n'.format(output))
    log.write('**********\n')
    log.write('Beginning calculation(s)...\n\n')
    i = 0
    df = pd.DataFrame()
    for selection in selections:
        try:
            print('Calculating RMSD for selection {}: {}\nStart time: {}\n'.format(i, selection,now()))
            log.write('Calculating RMSD for selection {}: {}\n'.format(i, selection))
            log.write('Start time: {}\n'.format(now()))
            sele = traj.top.select(selection)
            reference = traj[ref_index]
            rms = mdtraj.rmsd(traj, reference, atom_indices=sele, precentered=precentered)
            print('Completed calculation for selection{}: {}\nEnd time: {}\n'.format(i, selection,now()))
            log.write('Completed calculation for selection.\n')
            log.write('End time: {}\n'.format(now()))
            log.write('**********\n')
            df[selection] = rms
        except Exception as e:
            log.write('Error calculating RMSD for selection:\n{}\n'.format(e))
            log.write('Will attempt to calculate RMSD for more selections, if any.\n')
        i += 1
    if df.empty:
        log.write('No RMSD calculations were completed. Please check error(s) and inputs.\n')
        log.close()
        return df
    print('RMSD calculations complete. End time {}\n'.format(now()))
    log.write('RMSD calculations complete. End time {}\n'.format(now()))
    df.index = traj.time
    df.name = 'rmsd'
    df = write_xvg(output, df, 'rmsd')
    log.write('Wrote output to {}\n'.format(output))
    log.write('Job complete.')
    log.close()
    return (df, output)
    
def write_xvg(output, df, job):
    ptypes = {
        'rmsd':'timeseries'
    }
    f = open(output, 'w')
    f.write('# This file was created {}\n'.format(now()))
    f.write(f'# Created by: pymd.mdrun.{job}\n')
    f.write(f'@    title "{job.upper()}"\n')
    f.write(f'@    xaxis label "Time (ps)"\n')
    f.write(f'@    yaxis label "RMSD (nm)"\n')
    f.write(f'@TYPE xy\n')
    f.write(f'@PTYPE {ptypes[job]}\n')
    i = 0
    for column in df.columns:
        f.write('@ s{} legend "{}"\n'.format(i, column))
        i += 1
    f.close()
    with open(output, 'a') as f:
        df.to_csv(f, header=False, sep='\t')
    df = PostProcess.metadata(output, df=df)
    return df



