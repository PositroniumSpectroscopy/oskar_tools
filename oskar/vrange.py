#!python
""" vrange: find the vertical range of traces in hdf5 file

    Copyright (c) 2015, 2016 UNIVERSITY COLLEGE LONDON
    @author: Adam
"""
from __future__ import print_function, division
import sys
import os
import argparse
import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm
import oskar

def vrange_h5(data, name, **kwargs):
    """ find the vertical range for each squid. """
    progress = kwargs.get('progress', False)
    quiet = kwargs.get('quiet', False)
    verbose = kwargs.get('verbose', False)
    max_squid = kwargs.get('max_squid', -1)
    #initialise output DataFrame
    heads = np.array([name+"_range_reps", name+"_range_mean", name+"_range_std",
                      name+"_range_sem"])
    # get squid values
    squids = np.sort(np.array([int(key) for key in data.keys()]))
    if max_squid > 0:
        squids = squids[squids <= max_squid]
    avDF = pd.DataFrame(index=squids, columns=heads)
    avDF.index.name = 'SQUID'
    rtot = 0 # number of traces
    #average data
    if verbose:
        print('find vertical range')
    for sq in tqdm(squids, smoothing=0.1, unit=' squids', leave=True,
                   disable=bool(not progress)):
        squid = str(sq)
        if name in data[squid]:
            trace = np.array(data[squid][name])
            rep = np.shape(trace)[0]
            rtot = rtot + rep
            ranges = trace.max(axis=1)-trace.min(axis=1)
            avDF.loc[sq, name+"_range_reps"] = rep
            avDF.loc[sq, name+"_range_mean"] = ranges.mean()
            avDF.loc[sq, name+"_range_std"] = ranges.std()
            avDF.loc[sq, name+"_range_sem"] = ranges.std()/np.sqrt(rep)
    if verbose:
        print('remove empty columns')
    avDF = avDF.dropna(axis=1, how='all')
    if len(avDF.columns) > 0:
        if not quiet:
            tot = len(avDF.index)
            ntot = len(avDF.dropna().index)
            print(name, ":", "found the average vertical range of", str(rtot),
                  "traces, for", str(ntot), "out of", str(tot), "seqs.")
        return avDF
    else:
        raise IOError("no data found for " + name)

def vrange_loop_h5(data, name, vDF, uDF, **kwargs):
    """ find the vertical range of data grouped by var values (i.e., over loops). """
    progress = kwargs.get('progress', False)
    quiet = kwargs.get('quiet', False)
    verbose = kwargs.get('verbose', False)
    max_squid = kwargs.get('max_squid', -1)
    if verbose:
        print('find vertical range over loops')
    heads = np.array([name+"_range_reps", name+"_range_mean", name+"_range_std",
                      name+"_range_sem"])
    avDF = pd.DataFrame(index=uDF.index, columns=heads)
    avDF.index.name = 'VID'
    rtot = 0 # number of traces
    if verbose:
        print('import data for unique vars')
    for vid in tqdm(uDF.index, smoothing=0.1, unit=' vids',
                    leave=True, disable=bool(not progress)):
        squids = vDF[(vDF[uDF.columns] == uDF.loc[vid]).all(1)].index.values
        if max_squid > 0:
            squids = squids[squids <= max_squid]
        trace = oskar.h5_array(data, squids, name, ignore_missing=True)
        rep = np.shape(trace)[0]
        if rep > 0:
            rtot = rtot + rep
            ranges = trace.max(axis=1)-trace.min(axis=1)
            avDF.loc[vid, name+"_range_reps"] = rep
            avDF.loc[vid, name+"_range_mean"] = ranges.mean()
            avDF.loc[vid, name+"_range_std"] = ranges.std()
            avDF.loc[vid, name+"_range_sem"] = ranges.std()/np.sqrt(rep)
    if verbose:
        print('remove empty columns')
    avDF = avDF.dropna(axis=1, how='all')
    if len(avDF.columns) > 0:
        if not quiet:
            tot = len(avDF.index)
            ntot = len(avDF.dropna().index)
            ldat = len(data) if max_squid < 0 else np.min([max_squid, len(data)])
            nloops = float(ldat) / tot
            print(name, ":", "found the average vertical range of", str(rtot),
                  "traces, for", str(ntot), "out of", str(tot), "unique seqs from",
                  "%.2f" % nloops, "loop(s)")
        return avDF
    else:
        raise IOError("no data found for " + name)

def main():
    """ find the vertical range of traces in hdf5 file. """
    # defaults
    settings = oskar.Defaults()
    ddir = settings.values['dire'] if 'dire' in settings.values else None
    drid = settings.values['rid'] if 'rid' in settings.values else None
    dftype = settings.values['vrange'] if 'vrange' in settings.values else None
    # inputs
    # data options
    parser = argparse.ArgumentParser(description='Calculate the vertical range \
                            of waveforms (defaults.json key: vrange). Output \
                            the average and sem to *.dat (tsv) file(s).')
    dat_parser = parser.add_argument_group('HDF5 data')
    dat_parser.add_argument('-d', '--dire', nargs=1, default=ddir,
                            help='data directory, e.g. --dire \"Z:\\Data\"')
    dat_parser.add_argument('-r', '--rid', nargs=1, default=drid,
                            help='rid, e.g. --rid \"20160203_185633\"')
    dat_parser.add_argument('-f', '--ftype', nargs='+', default=dftype,
                            help='file type(s) to analyse, \
                            e.g. --ftype \"CH_L1"')
    dat_parser.add_argument('-L', '--loop', action="store_true", default=False,
                            help="average values over loops")
    dat_parser.add_argument('-M', '--max_squid', dest='max_squid', type=int, default=-1,
                            help="only analyse up to squid < max_squid")
    dat_parser.add_argument('-V', '--lvars', nargs='*', default=None,
                            help='specific vars to consider when looking for \
                            unique combinations (looping).')
    # script options
    script_parser = parser.add_argument_group('script options')
    script_parser.add_argument('-i', '--info', action="store_true", default=False,
                               help="pprint h5 file attributes")
    script_parser.add_argument('-b', '--progress', action="store_true", default=False,
                               help="display progress bar")
    script_parser.add_argument('-q', '--quiet', action="store_true", default=False,
                               help="surpress script output")
    script_parser.add_argument('-v', '--verbose', action="store_true", default=False,
                               help="verbose script information")
    # defaults.json
    def_parser = parser.add_argument_group('defaults.json')
    def_parser.add_argument('-s', '--set', action="store_true", default=False,
                            help="save args. (dire/ rid/ ftype) as default values")
    args = parser.parse_args()
    # data dire
    if args.dire is not None:
        # overwrite default dire
        settings.assign('dire', args.dire)
    else:
        raise Exception("data directory not specified, nor found in defaults.json. Use flag --dire")
    # rid
    if args.rid is not None:
        # overwrite default rid
        settings.assign('rid', args.rid)
    else:
        raise Exception("rid not specified, nor found in defaults.json. Use flag --rid")
    # ftype
    if args.ftype is not None:
        # overwrite default rid
        settings.assign('vrange', args.ftype)
        ftypes = args.ftype
    else:
        raise Exception("ftype not specified nor found in defaults.json (key: vrange). \
                         Use flag --ftype")
    if args.set:
        settings.save()
    # run averaging
    for rid in args.rid:
        if not (args.info or args.quiet):
            print(rid)
        if args.verbose:
            print('loading hdf5 file')
        h5 = oskar.H5Data(rid, args.dire[0])
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
                vDF = oskar.varDF(data)
                uDF = oskar.unqDF(vDF, lvars=args.lvars)
                if args.verbose:
                    print('save to unique_vars.dat')
                out_fil = os.path.join(out_dire, "unique_vars.dat")
                uDF.to_csv(out_fil, sep='\t')
            # average for each ftype
            for name in ftypes:
                sys.stdout.flush()
                if args.loop:
                    avDF = vrange_loop_h5(data, name, vDF, uDF, **vars(args))
                else:
                    avDF = vrange_h5(data, name, **vars(args))
                if args.verbose:
                    print('save average data')
                fil = "range_"+name+".dat"
                out_fil = os.path.join(out_dire, fil)
                if args.verbose:
                    print(out_fil)
                pd.DataFrame(avDF).to_csv(out_fil, sep='\t')
            if args.verbose:
                print("Done!")

if __name__ == "__main__":
    # execute only if run as a script
    main()
