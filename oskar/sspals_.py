#!python
""" sspals: Combine hdf5 data (chmx).  Re-analyse each to find t0 (cfd
    trigger) and the delayed fraction (DF = BC/ AC) for limits (A,B,C).

    Copyright (c) 2015, 2016 UNIVERSITY COLLEGE LONDON
    @author: Adam Deller
"""
from __future__ import print_function, division
import sys
import os
import re
import argparse
import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm
import sspals
import oskar

def sspals_h5(data, hi, low, **kwargs):
    """ sspals analysis for each squid. """
    progress = kwargs.get('progress', False)
    quiet = kwargs.get('quiet', False)
    verbose = kwargs.get('verbose', False)
    max_squid = kwargs.get('max_squid', -1)
    #initialise output DataFrame
    columns = ['t0', 'DF']
    heads = np.array([[k+"_reps", k+"_mean", k+"_std",
                       k+"_sem"] for k in columns]).flatten()
    # get squid values
    squids = np.sort(np.array([int(key) for key in data.keys()]))
    if max_squid > 0:
        squids = squids[squids <= max_squid]
    av_df = pd.DataFrame(index=squids, columns=heads)
    av_df.index.name = 'SQUID'
    rtot = 0 # number of repeats
    #average data
    if verbose:
        print('build chmx from', hi, 'and', low, 'data, and analyse.')
    for sq in tqdm(squids, smoothing=0.1, unit=' squids', leave=True,
                   disable=bool(not progress)):
        squid = str(sq)
        if (hi in data[squid]) and (low in data[squid]):
            hdat = np.array(data[squid][hi])
            ldat = np.array(data[squid][low])
            dt = oskar.arr2dict(data[squid][hi].attrs.items())['dt']
            chmx = sspals.chmx(hdat, ldat, **kwargs)
            ntraces = np.shape(chmx)[0]
            if ntraces > 0:
                rtot = rtot + ntraces
                sspalsDF = pd.DataFrame(sspals.sspals(chmx, dt, dropna=True, **kwargs))
                for col in columns:
                    if col in sspalsDF:
                        rep = sspalsDF[col].count()
                        av_df.loc[sq, col+'_ra_reps'] = rep
                        av_df.loc[sq, col+'_ra_mean'] = sspalsDF[col].mean()
                        av_df.loc[sq, col+'_ra_std'] = sspalsDF[col].std()
                        av_df.loc[sq, col+'_ra_sem'] = sspalsDF[col].std()/np.sqrt(rep)
    if verbose:
        print('remove empty columns')
    av_df = av_df.dropna(axis=1, how='all')
    if len(av_df.columns) > 0:
        if not quiet:
            tot = len(av_df.index)
            ntot = len(av_df.dropna().index)
            print("SSPALS [", hi, ', ', low, "] : re-analysed ",
                  str(rtot), " spectra, for", str(ntot), "out of", str(tot), "sequences.")
        return av_df
    else:
        raise IOError("no data found for " + np.array_str([hi, low]))

def sspals_loop_h5(data, hi, low, v_df, u_df, **kwargs):
    """ sspals analysis for data grouped by var values (i.e., over loops). """
    progress = kwargs.get('progress', False)
    quiet = kwargs.get('quiet', False)
    verbose = kwargs.get('verbose', False)
    max_squid = kwargs.get('max_squid', -1)
    if verbose:
        print('find vertical range over loops')
    columns = ['t0', 'DF']
    heads = np.array([[k+"_reps", k+"_mean", k+"_std",
                       k+"_sem"] for k in columns]).flatten()
    av_df = pd.DataFrame(index=u_df.index, columns=heads)
    av_df.index.name = 'VID'
    rtot = 0 # number of traces
    if verbose:
        print('import data for unique vars')
    for vid in tqdm(u_df.index, smoothing=0.1, unit=' vids',
                    leave=True, disable=bool(not progress)):
        squids = v_df[(v_df[u_df.columns] == u_df.loc[vid]).all(1)].index.values
        if max_squid > 0:
            squids = squids[squids <= max_squid]
        hdat = oskar.h5_array(data, squids, hi, ignore_missing=True)
        ldat = oskar.h5_array(data, squids, low, ignore_missing=True)
        dt = oskar.arr2dict(data[str(squids[0])][hi].attrs.items())['dt']
        chmx = sspals.chmx(hdat, ldat, **kwargs)
        ntraces = np.shape(chmx)[0]
        if ntraces > 0:
            rtot = rtot + ntraces
            sspalsDF = pd.DataFrame(sspals.sspals(chmx, dt, dropna=True, **kwargs))
            for col in columns:
                if col in sspalsDF:
                    rep = sspalsDF[col].count()
                    av_df.loc[vid, col+'_ra_reps'] = rep
                    av_df.loc[vid, col+'_ra_mean'] = sspalsDF[col].mean()
                    av_df.loc[vid, col+'_ra_std'] = sspalsDF[col].std()
                    av_df.loc[vid, col+'_ra_sem'] = sspalsDF[col].std()/np.sqrt(rep)
    if verbose:
        print('remove empty columns')
    av_df = av_df.dropna(axis=1, how='all')
    if len(av_df.columns) > 0:
        if not quiet:
            tot = len(av_df.index)
            ntot = len(av_df.dropna().index)
            ldat = len(data) if max_squid < 0 else np.min([max_squid, len(data)])
            nloops = float(ldat) / tot
            print("SSPALS [", hi, ', ', low, "] : re-analysed ",
                  str(rtot), " spectra, for", str(ntot), "out of", str(tot), "unique seqs from",
                  "%.2f" % nloops, "loop(s)")
        return av_df
    else:
        raise IOError("no data found for " + np.array_str([hi, low]))

def main():
    """ sspals analysis of traces in hdf5 file. """
    # defaults
    settings = oskar.Defaults()
    bdir = settings.values['base'] if 'base' in settings.values else None
    ddir = settings.values['dire'] if 'dire' in settings.values else None
    drid = settings.values['rid'] if 'rid' in settings.values else None
    dftype = settings.values['sspals'] if 'sspals' in settings.values else None
    # inputs
    # data options
    parser = argparse.ArgumentParser(description='Calculate the delayed fraction \
                            of waveforms made by combining two (defaults.json \
                            key: sspals), applying a constant-fraction discriminator, \
                            and integrating. Output the average and sem to *.dat (tsv) file(s).')
    dat_parser = parser.add_argument_group('HDF5 data')
    dat_parser.add_argument('-b', '--base', nargs=1, default=bdir,
                            help='base directory, e.g. --base \"Z:\\Data\"')
    dat_parser.add_argument('-d', '--dire', nargs=1, default=ddir,
                            help='data directory. Defaults to "[base]\\YYYY\\mm\\dd\\rid"')
    dat_parser.add_argument('-r', '--rid', nargs=1, default=drid,
                            help='rid, e.g. --rid \"20160203_185633\"')
    dat_parser.add_argument('-f', '--ftype', nargs=2, default=dftype,
                            help='file types [hi, low] to combine into chmx data, \
                            e.g. --ftype "CH_L0" "CH_L1"')
    dat_parser.add_argument('-L', '--loop', action="store_true", default=False,
                            help="average values over loops")
    dat_parser.add_argument('-M', '--max_squid', dest='max_squid', type=int, default=-1,
                            help="only analyse up to squid < max_squid")
    dat_parser.add_argument('-V', '--lvars', nargs='*', default=None,
                            help='specific vars to consider when looking for \
                            unique combinations (looping).')

    # sspals options
    sspals_parser = parser.add_argument_group('SSPALS settings')
    sspals_parser.add_argument('--validate', action="store_true", default=False,
                               help="validate data.")
    sspals_parser.add_argument('--min_range', dest='min_range', type=float, default=0.1,
                               help="minimum vertical range for validation, default 0.1")
    sspals_parser.add_argument('--n_bsub', dest='n_bsub', type=int, default=100,
                               help="number of data points to use for background subtraction,\
                               default 100")
    sspals_parser.add_argument('--limits', dest='limits', type=str, default="-1e-8, 3.5e-8, 6.5e-7",
                               help="integration bounds, --limits \"A B C\".")
    # cfd options
    cfd_parser = parser.add_argument_group('CFD settings')
    cfd_parser.add_argument('--cfd_scale', dest='cfd_scale', type=float, default=0.8,
                            help="cfd scale, default 0.8")
    cfd_parser.add_argument('--cfd_offset', dest='cfd_offset', type=float, default=1.4E-8,
                            help="cfd offset, default 1.4E-8")
    cfd_parser.add_argument('--cfd_threshold', dest='cfd_threshold', type=float, default=0.04,
                            help="cfd threshold, default 0.04")
    # script options
    script_parser = parser.add_argument_group('script options')
    script_parser.add_argument('-i', '--info', action="store_true", default=False,
                               help="pprint h5 file attributes")
    script_parser.add_argument('-t', '--progress', action="store_true", default=False,
                               help="display tqdm progress bar")
    script_parser.add_argument('-q', '--quiet', action="store_true", default=False,
                               help="surpress script output")
    script_parser.add_argument('-v', '--verbose', action="store_true", default=False,
                               help="verbose script information")
    # defaults.json
    def_parser = parser.add_argument_group('defaults.json')
    def_parser.add_argument('-s', '--set', action="store_true", default=False,
                            help="save args. (dire/ rid/ ftype) as default values")
    args = parser.parse_args()
    # data
    if args.base is not None:
        # overwrite default base directory
        settings.assign('base', args.base)
    else:
        raise Exception("base directory not specified, nor found in defaults.json. Use flag --base")
    if args.dire is not None:
        # overwrite default data directory
        settings.assign('dire', args.dire)
    else:
        args.dire = [None]
    # rid
    if args.rid is not None:
        # overwrite default rid
        settings.assign('rid', args.rid)
    else:
        raise Exception("rid not specified, nor found in defaults.json. Use flag --rid")
    # ftype
    if args.ftype is not None:
        # overwrite default rid
        settings.assign('sspals', args.ftype)
        ftypes = args.ftype
        [hi, low] = ftypes
    else:
        raise Exception("ftype not specified nor found in defaults.json (key: vrange). \
                         Use flag --ftype")
    if args.set:
        settings.save()
    args.limits = np.array(re.split(",", args.limits), dtype=float)
    # run averaging
    for rid in args.rid:
        if not (args.info or args.quiet):
            print(rid)
        if args.verbose:
            print('loading hdf5 file')
        h5 = oskar.H5Data(rid, args.base[0], args.dire[0])
        if args.info:
            h5.pprint()
        # output
        if args.verbose:
            print('building output path')
        out_folder = 'Average_Loops' if args.loop else 'Average'
        out_dire = h5.out_dire(out_folder)
        if args.verbose:
            print(out_dire)
        # load hdf5 file
        with h5py.File(h5.fil, 'r') as dfil:
            data = dfil['.']
            if args.loop:
                if args.verbose:
                    print('finding unique var combinations')
                # read vars from metadata
                v_df = oskar.var(data)
                u_df = oskar.unique(v_df, lvars=args.lvars)
                if args.verbose:
                    print('save to unique_vars.dat')
                out_fil = os.path.join(out_dire, "unique_vars.dat")
                u_df.to_csv(out_fil, sep='\t')
            sys.stdout.flush()
            if args.loop:
                av_df = sspals_loop_h5(data, hi, low, v_df, u_df, **vars(args))
            else:
                av_df = sspals_h5(data, hi, low, **vars(args))
            if args.verbose:
                print('save average data')
            fil = "sspals_"+hi+ np.array_str(args.limits, precision=2) + ".dat"
            out_fil = os.path.join(out_dire, fil)
            if args.verbose:
                print(out_fil)
            pd.DataFrame(av_df).to_csv(out_fil, sep='\t')
            if args.verbose:
                print("Done!")

if __name__ == "__main__":
    # execute only if run as a script
    main()
