import marimo

__generated_with = "0.11.31"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # scATAC quality control using snapATAC2
        This notebook will guide you through basic quality control of scATAC data using snapATAC2.
        """
    )
    return


@app.cell(hide_code=True)
def imports():
    from dask.diagnostics import ProgressBar
    from pathlib import Path

    import io
    import os
    import json
    import logging
    import rapidgzip
    import dask.dataframe as dd
    import numpy as np
    import pandas as pd
    import pyarrow as pa

    # Initialize logger.
    logger = logging.getLogger(__name__)
    # Configure logging with a timestamp
    logging.basicConfig(
        filename='logs/marimo_snapatac2.log.txt',
        level=logging.DEBUG,
        encoding='utf-8',
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"  # Human-readable timestamp format
    )

    os.environ["SSL_CERT_FILE"] = "/Users/emattei/GitHub/epi-SHARE-seq-pipeline/.venv/lib/python3.10/site-packages/certifi/cacert.pem"

    # Enable the progress bar
    pbar = ProgressBar()
    pbar.register()
    return (
        Path,
        ProgressBar,
        dd,
        io,
        json,
        logger,
        logging,
        np,
        os,
        pa,
        pbar,
        pd,
        rapidgzip,
    )


@app.cell(hide_code=True)
def _(mo):
    analysis_parameters = mo.query_params()
    return (analysis_parameters,)


@app.cell
def parameters(analysis_parameters, json, logging):
    logging.info(json.dumps(str(analysis_parameters), indent=4))
    analysis_parameters
    return


@app.cell(hide_code=True)
def fragment_form_definition(Path, mo):
    # Function to validate to the fragment file path to make sure it exists and is a file.
    def validate_fragment_form(input):
        return None if Path(input).exists() else "ERROR: File does not exist."

    fragment_file_form=mo.ui.text(
                                  full_width=True,
                                ).form(
                                    bordered=False,
                                    submit_button_label="Load fragments",
                                    validate=validate_fragment_form,
                                )
    return fragment_file_form, validate_fragment_form


@app.cell(hide_code=True)
def fragment_form_input(fragment_file_form, mo):
    mo.vstack([
        mo.md(
            '''
            **::lucide:dna:: Enter fragment file path:**
            '''
        ),
        fragment_file_form
    ])
    return


@app.cell
def _():
    path = "src\\python\\qc_atac\\data\\IGVFDS4145LQGK.fragments.tsv.bgz"
    return (path,)


@app.cell
def _(io, path, rapidgzip):
    bc_count = set()
    region_count = set()
    fragment_count = 0
    with rapidgzip.open(path, parallelization=8) as fh:
        with io.TextIOWrapper(fh) as f:
            for line in f:
                col = line.strip().split("\t")
                bc_count.add(col[4])
                region_count.add(f"{col[0]}_{str(col[1])}_{str(col[2])}")
                fragment_count += 1

    len(bc_count)
    len(region_count)
    fragment_count
    return bc_count, col, f, fh, fragment_count, line, region_count


@app.cell
def dask_ingestion(dd, np, path):
    blocksize = 64 * 1024 * 1024  # 64 MB block size

    temp_dask_df = dd.read_csv(path,
                               sep="\t",
                               names=["chr", "start", "end", "barcode", "supporting_reads"],
                               dtype={"chr": "category",
                                      "start": np.uint32,
                                      "end": np.uint32,
                                      "barcode": "string",
                                      "supporting_reads": np.uint16,
                                     },
                               compression="gzip"
                              )
    return blocksize, temp_dask_df


@app.cell
def _(temp_dask_df):
    print(type(temp_dask_df))  # Should print <class 'dask.dataframe.core.DataFrame'>
         # Should print <class 'dask.dataframe.core.DataFrame'>
    return


@app.cell
def _(temp_dask_df):
    import zarr
    temp_dask_df.to_dask_array(lengths=True).to_zarr("src\\python\\qc_atac\\data\\IGVFDS4145LQGK",
                        write_empty_chunks=False,
                            mode="w",
                           )

    final_path = "src\python\qc_atac\data\IGVFDS4145LQGK"
    return final_path, zarr


@app.cell
def _(dd, final_path):
    ddf = dd.read_zarr(final_path)
    return (ddf,)


@app.cell
def raw_metrics(ddf, format_number, logging, mo):
    mo.stop(True)
    logging.info(f"Computing statistics.")
    # Get the total number or rows in the fragment file.
    total_rows = ddf.shape[0].compute()
    # Get the number of unique barcodes.
    unique_barcodes = ddf["barcode"].nunique().compute()
    # Get the total number of molecules(unique+duplicates).
    total_molecules = ddf["supporting_reads"].sum().compute()

    number_of_regions = mo.stat(
        value=f"{format_number(total_rows)}",
        label="Number of fragments",
    )

    number_of_barcodes = mo.stat(
        value=f"{format_number(unique_barcodes)}",
        label="Number of unique barcodes",
    )

    percent_duplicates = mo.stat(
        value=f"{(total_molecules-total_rows)/total_molecules:.1%}", 
        label="Percent duplicates",
    )
    logging.info(f"Computing statistics.-DONE")
    mo.hstack([number_of_regions, number_of_barcodes, percent_duplicates], justify="center", gap="2rem")
    return (
        number_of_barcodes,
        number_of_regions,
        percent_duplicates,
        total_molecules,
        total_rows,
        unique_barcodes,
    )


@app.cell
def _(mo):
    mo.md(
        r"""
        ## Check for over-represented genomic regions across barcodes.
        We want to check if there are regions that are associated with an unusually high number of barcodes.
        This can indicate potential issues with the library preparation, mappability issues, or PCR amplification artifacts.
        """
    )
    return


@app.cell
def _():
    #barcode_counts_per_region = ddf.groupby(["chr", "start", "end"])["barcode"].nunique(split_out=12).compute()
    return


@app.cell
def _(ddf, mo, pd):
    mo.stop(True)
    def compute_metrics():
        # Each line is chr, start, end, bracode, duplicates
        print("Computing metrics")
        return pd.DataFrame({
            "chrom": ["chr1", "chr2", "chr3"],
            "start": [100, 200, 300],
            "end": [200, 300, 400],
            "barcode_count": [100, 200, 300],
        })
    # Trigger computation (e.g., count rows)
    row_count = ddf.shape[0].compute()
    return compute_metrics, row_count


@app.cell
def _(ddf, mo):
    mo.stop(True)
    # Filter by barcodes of interest
    barcodes_of_interest = {"GGAATGATGACGGATTGTGTCCTT_IGVFSM4419MRHA", "TACCGAGCACTTGATGCGATGTTT_IGVFSM4419MRHA"}
    filtered = ddf[ddf["barcode"].isin(barcodes_of_interest)]

    # Compute the result (triggers execution)
    result = filtered.compute()
    return barcodes_of_interest, filtered, result


@app.cell
def format_number():
    def format_number(num):
        if num >= 1_000_000_000:
            return f"{num / 1_000_000_000:.1f}B"  # Billions
        elif num >= 1_000_000:
            return f"{num / 1_000_000:.1f}M"  # Millions
        elif num >= 1_000:
            return f"{num / 1_000:.1f}K"  # Thousands
        else:
            return str(num)  # Less than 1,000
    return (format_number,)


@app.cell
def _(Path):
    def get_parquet_file_path(input_path):
        try:
            prefix = Path(input_path)
        except Exception as e:
            return f"ERROR: {e}"
        while prefix.suffix:
            prefix = prefix.with_suffix("")
        prefix = prefix.stem
        return Path(input_path).parent / f"{prefix}.parquet"
    return (get_parquet_file_path,)


if __name__ == "__main__":
    app.run()
