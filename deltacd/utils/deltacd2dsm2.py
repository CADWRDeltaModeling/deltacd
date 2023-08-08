# -*- coding: utf-8 -*-
"""
Postprocessing tools for deltacd
"""
from pathlib import Path
import pandas as pd
import xarray as xr
import pyhecdss


def deltacd2dsm2(ncfn, dssfn, inpfn=None):
    """
    Converting deltacd output from nc to dss format

    Parameters
    ----------
    ncfn : str
        netcdf filename.nc from DCD
    dssfn : str
        output dss filename.dss
    inpfn : str
        output dsm2 sources and sinks input filename.inp

    Returns
    -------
    None.

    """
    deltacd = xr.open_dataset(str(ncfn))

    dcd = pyhecdss.DSSFile(str(dssfn), create_new=True)
    ptype = "PER-AVER"
    units = "CFS"
    A = "DELTACD-HIST+NODE"
    E = "1DAY"
    F = "DWR-BDO"

    if inpfn is not None:
        # creating corresponding source flow input
        Descriptions = "#Delta Island Consumptive Use Estimated by Delta CD"
        header = "\nSOURCE_FLOW"
        with open(inpfn, "w") as writer:
            writer.writelines(Descriptions)
            writer.writelines(header)
            writer.writelines("\nNAME           NODE SIGN FILLIN FILE        PATH")
        FILLIN = "last"
        FILE = "${DICUFILE}"

    res_lines = []
    for n in deltacd.node.values:
        for v in list(deltacd.keys()):
            B = n
            C = v
            outpath = "/%s/%s/%s//%s/%s/" % (A, B, C, E, F)
            ts = deltacd[v].sel(node=n).to_series()
            if (ts == 0).all():  # ignoring the node that has zero flow values.
                print("node %s has zero flow: %s" % (n, v))
                # continue
            if ts.index.freq is None:
                ts.index.freq = ts.index.inferred_freq
            dcd.write_rts(outpath, ts, units, ptype)
            print(outpath + "written!")

            if inpfn is not None:
                NAME = "dicu_%s_%s" % (v.split("-")[0].lower(), n)
                NAME = NAME.ljust(16)
                NODE = n.rjust(3)
                SIGN = -1
                if v == "DRAIN-FLOW":
                    SIGN = 1
                if n == "BBID":  # this is reservoir
                    newline = "\n%s%s%5s%5s%14s %s" % (
                        NAME,
                        "clifton_court",
                        SIGN,
                        FILLIN,
                        FILE,
                        outpath,
                    )
                    res_lines.append(newline)
                else:
                    newline = "\n%s%s%5s%5s%14s %s" % (
                        NAME,
                        NODE,
                        SIGN,
                        FILLIN,
                        FILE,
                        outpath,
                    )
                    with open(inpfn, "a") as writer:
                        writer.writelines(newline)
    dcd.close()
    deltacd.close()
    print("dss file written!")

    if inpfn is not None:
        with open(inpfn, "a") as writer:
            writer.writelines(
                [
                    "\nEND",
                    "\n",
                    "\nSOURCE_FLOW_RESERVOIR",
                    "\nNAME            RES_NAME      SIGN FILLIN FILE        PATH",
                ]
            )
            writer.writelines(res_lines)
            writer.writelines("\nEND")
        print("dsm2 input file written!")

    return


def time_index_correction(ts_series):
    time_str = [str(d) for d in ts_series.index]
    time_index = pd.to_datetime(time_str)
    modified_ts = pd.Series(ts_series.values, index=time_index).resample("1D").nearest()
    return modified_ts


def read_dcd(dcd_fn, nodes, start_time, end_time):
    """
    Parameters
    ----------
    dcd_fn : str
        DCD input dss file for dsm2.
    nodes : list
        A list of node numbers.
    start_time : str
        eg., date in the format of "2020-01-12"
    end_time : str
        Similar as above

    Returns
    -------
    drain_df : pandas dataseries
        drainage
    div_df: pandas dataseries
        diversion
    seep_df: pandas dataseries
        seepage
    """

    d = pyhecdss.DSSFile(dcd_fn)
    cat = d.read_catalog()
    nodes_cat = cat[cat.B.isin(nodes)]
    cat_drain = nodes_cat[nodes_cat.C.isin(["DRAIN-FLOW"])]
    cat_div = nodes_cat[nodes_cat.C.isin(["DIV-FLOW"])]
    cat_seep = nodes_cat[nodes_cat.C.isin(["SEEP-FLOW"])]
    # cat_div = nodes_cat[nodes_cat.C.isin(['DIV-FLOW','SEEP-FLOW'])]

    # calculate the total flow of these nodes
    paths_drain = d.get_pathnames(cat_drain)
    df = pd.DataFrame()
    for b, c, p in zip(cat_drain.B, cat_drain.C, paths_drain):
        data = d.read_rts(p, start_time, end_time)
        ts = data.data[data.data.columns[0]]
        df["%s_%s" % (b, c)] = ts

    df_drain = df.sum(axis=1)
    df_drain = time_index_correction(df_drain)

    paths_div = d.get_pathnames(cat_div)
    df = pd.DataFrame()
    for b, c, p in zip(cat_div.B, cat_div.C, paths_div):
        data = d.read_rts(p, start_time, end_time)
        ts = data.data[data.data.columns[0]]
        df["%s_%s" % (b, c)] = ts

    df_div = df.sum(axis=1)
    df_div = time_index_correction(df_div)

    paths_seep = d.get_pathnames(cat_seep)
    df = pd.DataFrame()
    for b, c, p in zip(cat_seep.B, cat_seep.C, paths_seep):
        data = d.read_rts(p, start_time, end_time)
        ts = data.data[data.data.columns[0]]
        df["%s_%s" % (b, c)] = ts

    df_seep = df.sum(axis=1)
    df_seep = time_index_correction(df_seep)
    d.close()

    return df_drain, df_div, df_seep


def read_net_div(dcd_fn, nodes, start_time, end_time):
    """
    Parameters
    ----------
    dcd_fn : str
        DCD input dss file for dsm2.
    nodes : list
        A list of node numbers.
    start_time : str
        eg., date in the format of "2020-01-12"
    end_time : str
        Similar as above

    Returns
    -------
    net_df : pandas series
        Total diversion rate

    """
    drain, div, seep = read_dcd(dcd_fn, nodes, start_time, end_time)
    net_df = div - drain + seep
    return net_df


def create_argparser():
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert DeltaCD output to DSM2 DSS input files"
    )
    parser.add_argument(
        "--input",
        dest="input",
        type=Path,
        required=True,
        help="DeltaCD output netcdf filename (input file)",
    )
    parser.add_argument(
        "--output_dss",
        dest="output_dss",
        type=Path,
        required=True,
        help="DSM2 DSS filename (output file)",
    )
    parser.add_argument(
        "--output_inp",
        dest="output_inp",
        type=Path,
        help="DCD *.inp filename for DSM2 (output file)",
    )
    args = parser.parse_args()
    return args


def main():
    args = create_argparser()
    deltacd2dsm2(args.input, args.output_dss, args.output_inp)


if __name__ == "__main__":
    main()
