# TODO: manual test:
# /home/mvn373/fast/tools/python36/bin/python3 -u /home/mvn373/fast/prs-workshop/apply-pgscatalog-prs.py --genetic /home/mvn373/fast/prs-repo/pre-calculated/genetic/TARGET_final --prs-wm /home/mvn373/fast/prs-repo/pre-calculated/WMs/olink/PGS000217.txt --out /home/mvn373/fast/prs-repo/pre-calculated/scores/tool/dev-test
# diff /home/mvn373/fast/prs-repo/pre-calculated/scores/tool/olink000217-TARGETold.prs /home/mvn373/fast/prs-repo/pre-calculated/scores/tool/dev-test.prs
# TODO: make support of ORs
import argparse
from io import TextIOWrapper
import logging
import time
import os
import shutil
import subprocess
from typing import Union

import pandas as pd
import yaml

from rich.align import Align
from rich.console import Console
from rich.layout import Layout
from rich.live import Live
from rich.panel import Panel
from rich.prompt import Confirm, Prompt
from rich.table import Table

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--bfile", "--genetic", default=None, type=str,
                    help="Absolute path prefix to plink19 trio bed/bim/fam files "
                         "(for example, for the file '/emc/data/File.bed' this would be '/emc/data/File')")
parser.add_argument("--pfile", default=None, type=str,
                    help="Absolute path prefix to plink2 trio pgen/pvar/psam files "
                         "(for example, for the file '/emc/data/File.pgen' this would be '/emc/data/File')")
parser.add_argument("-p", "--prs-wm", default=None, type=str,
                    help="Absolute path to PRS weight matrix downloaded from PGSCatalog")
parser.add_argument("-o", "--out", default=None, type=str,
                    help="Output path prefix "
                         "(for example, to generate the files '/home/abc123/prs/output.profile' "
                         "this would be '/home/abc123/prs/output')")
parser.add_argument("-c", "--config", default=None, type=str,
                    help="Path to yaml config")
arguments = parser.parse_args()

console = Console()
layout = Layout()
live = Live(console=console)
CONFIG: Union[None, dict] = None
LOGFILE: Union[None, TextIOWrapper] = None


PLINK_KEY_COLUMN = "rsid"
PGSCATALOG_KEY_COLUMN = "pos_key"


def check_input(args):
    global PLINK_FILEFORMAT
    if args.bfile is not None:
        args.genetic = args.bfile
        PLINK_FILEFORMAT = ("bed", "bim", "fam")
    elif args.pfile is not None:
        args.genetic = args.pfile
        PLINK_FILEFORMAT = ("pgen", "pvar", "psam")
    else:
        print(f"plink genetic data are required. Please specify either a plink19 or plink2 dataset.")
        exit(4)
    for ext in PLINK_FILEFORMAT:
        if not os.path.isfile(f"{args.genetic}.{ext}"):
            print(f"plink genetic data are malformed - {args.genetic}.{ext} is missing")
            exit(4)
    if not hasattr(args, 'genetic') or args.genetic is None:
        print("plink genetic data are required")
        parser.print_help()
        exit(1)
    if args.prs_wm is None:
        print("PRS weight matrix is required")
        parser.print_help()
        exit(2)
    if not os.path.isfile(args.prs_wm):
        print(f"PRS weight matrix file does not exist")
        exit(5)
    if args.out is None:
        print("Output prefix is required")
        parser.print_help()
        exit(3)


def maybe_load_config(config_file: Union[str, None]):
    global CONFIG

    if config_file is None:
        return
    with open(config_file, "r") as yamlfile:
        CONFIG = yaml.load(yamlfile, Loader=yaml.FullLoader)
        printout(f"===============================")
        printout(f"CONFIG WAS PROVIDED AT {config_file}")
        printout(f"RUNNING IN AUTOMATED MODE")
        printout(f"===============================")


def ask(text: str, _tool, _config: str, **kwargs):
    if CONFIG is None:
        # Manual mode, ask user
        return _tool.ask(text, **kwargs)
    else:
        # Automated mode
        if _config not in CONFIG:
            raise RuntimeError(f"An answer to the question {text} should have been provided in the config as {_config}, "
                               f"but there was no value for this parameter")
        return CONFIG[_config]


def print_error_files():
    is_to_print_the_guide = ask("Would you like to see a troubleshoot guide?", _config="n", default="n", _tool=Confirm)
    if is_to_print_the_guide:
        print("ERROR: 1A. Either the plink data don't use the conventional up-to-date rsID.")
        print("ERROR: 1B. Or the PGSCatalog data don't use the conventional up-to-date rsID.")
        print("ERROR: 2A. Or you are using un-imputed plink data. THEN you have to impute the data.")
        print("ERROR: *2B. Or you have data imputed to a very specific panel that does not contain a lot of common variants"
              " (happens seldom, probably won't happen to you). THEN you have to impute the data to a more general panel.")
        print("ERROR: 3. Or, if you are using chromosome:position fallback, the PRS weight matrix might be using "
              "a different reference genome than your plink data. "
              "THEN you need to liftover your data to the correct reference genome.")
        print("ERROR: 4. Or the PRS was done using a very specific chip and contains such unique SNPs "
              "that no imputation panel has them (happens seldom, probably won't happen to you). "
              "THEN you have to find an imputation panel that has these SNPs and impute the data to this panel.")
        print("Here is your troubleshoot guide:")
        print("---")
        print("1A")
        print("Q: How do I figure out if plink data are using proper rsID? (problem 1A)")
        print("A: Ask your data provider (typically the phenomics platform) about it. "
              "Or, look at the files - you might notice non-'rsNNNNNN' formatted ids (a few dots are ok).")
        print("Q: What can I do about plink data rsID annotations?")
        print("A: Ask your data provider (typically the phenomics platform) to annotate the rsID in your data.")
        print("---")
        print("1B")
        print("Q: How do I figure out if PGSCatalog data are using proper rsID? (problem 1B)")
        print("A: Look at the files - you might notice non-'rsNNNNNN' formatted ids (a few dots are ok).")
        print("Q: What can I do about PGSCatalog data rsID annotations?")
        print("A: Try using a fallback 'chromosome:position' method of matching. "
              "Remember, you will need to replace the ID in plink 'bim' file.")  # TODO: make it an option in the code
        print("---")
        print("2A, *2B")
        print("Q: How do I figure out if the data are imputed or not, and to which panel? (problems 2A, *2B)")
        print("A: Ask your data provider (typically the phenomics platform) about it. "
              "Imputed data will typically contain at least a few millions SNPs.")
        print("---")
        print("3")
        print("Q: How do I figure out if the reference genomes are different between the data? (problem 3)")
        print("A: Ask your data provider (typically the phenomics platform) which reference genome was used for your data. "
              "Check the reference genome of your score in the 'PGS Catalog Metadata' Excel sheet "
              "at the http://www.pgscatalog.org/downloads/ . See, if they match or not.")
        print("---")
        print("4")
        print("Q: How do I now if problem 4 is the case?")
        print("A: Unless you have a clear proof, it is almost definitely not.")
        print("---")
        print("General questions")
        print("Q: I have no idea about reference genomes and imputation, what do I do?")
        print("A: We did our best we could do automatically processing the PGSCatalog data. "
              "But if want to run the analysis, it is ultimately your responsibility to assure that you have proper data. "
              "You don't want to publish the results that are fake because you mishandled the data, do you? ")
        print("Ideally, your data provider (typically the phenomics platform) would *always* supply data to you "
              "in the same reference genome, imputed to the same imputation panel. If this was the case, this script would "
              "have been smarter at handling these issues automatically. But at the moment when this script was written "
              "this was not the case =/ Unfortunately, this means it is your problem now to figure it all out "
              "or to force your data provider to work properly ¯\_(ツ)_/¯")
        print("Q: I still have no idea, but I need this PRS, what do I do?")
        print("A: Invite a bioinformatician with programming / data science background, explain your problem to them, "
              "and show this error message.")


def print_non_annotated_plink_file_error():
    printout("ERROR: Your plink data have no variants IDs at all.")
    printout("ERROR: Ask your data supplier (probably the phenomics platform) to annotate it using dbSNP,")
    printout("ERROR: or at least assign the dummy 'chrom:pos' ids to the SNPs, if you know what you are doing.")
    printout("ERROR: This script do *not* need *dbSNP* variant IDs, but it *does* need *some* IDs.")
    printout("---")
    printout("P.S. We do *not* recommend using any kind of non-dbSNP ids,")
    printout(" because your colleagues' scripts may, unfortunately, rely on dbSNP annotation.")
    printout(" If you are ok with having 'chrom:pos' encoding and you think that you know what you are doing,")
    printout(" then you must be able to construct the annotation command using bash 'cut' and 'paste' yourself.")


# def input_files_qc(pgscatalog_df: pd.DataFrame, plink_variants_df: pd.DataFrame):
#     # Remove indels - we are unsure if they are dangerous, but we will remove them just in case
#     lines_were = pgscatalog_df.shape[0]
#     nucleotides = ("A", "C", "G", "T")
#     df = pgscatalog_df[pgscatalog_df["effect_allele"].isin(nucleotides)]
#     if "reference_allele" in df.columns:
#         df = df[df["reference_allele"].isin(nucleotides)]
#     printout(f"  {df.shape[0]} SNPs were preserved, {lines_were - df.shape[0]} indels were filtered")
#
#     # Chromosome names cutting - changes dataframes inplace
#     chrom_cutter = lambda x: x[3:] if x.startswith("chr") else x
#     pgscatalog_df["chr_name"] = pgscatalog_df["chr_name"].astype(str).apply(chrom_cutter)
#     plink_variants_df["chr_name"] = plink_variants_df["chr_name"].astype(str).apply(chrom_cutter)


GENERAL_DATA_FORMAT = """
pgscatalog_df: A dataframe of PGSCatalog score weight matrix data with columns:
            | chr_name | chr_position | effect_allele | reference_allele [may be omitted] | effect_weight |
            |----------|--------------|---------------|-----------------------------------|---------------|
            |     1    |   5894332    |       T       |                 G                 |   -0.000073   |
                    
plink_variants_df: A dataframe of local plink variants data with columns:
            | chr_name | rsid | chr_position | base_allele [A1] | alternative_allele [A2] |
            |----------|------|--------------|------------------|-------------------------|
            |     1    | rs12 |   5894332    |         T        |            G            |
"""


def check_plink_exec():
    # Check config for plink executable, else ask
    global PLINK_EXEC
    global PLINK_VERSION
    default = ["plink2"] + subprocess.run([
        "which" , "plink", "plink2"
    ], universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.rstrip().split("\n")
    PLINK_EXEC = ask(
        f"What path to the [bold red]plink executable[/bold red] should be used?",
        default=default[-1] if default[-1] else "plink2",
        _config="plink executable",
        _tool=Prompt,
    )
    PLINK_VERSION = int(subprocess.run([
        PLINK_EXEC, "--version"
    ], universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout[7])


def prep_keys_rsid(pgscatalog_df: pd.DataFrame, plink_variants_df: pd.DataFrame):
    # QC the plink rsid. Are rsids even present, or all rsids are just dots (".")?
    if plink_variants_df.head(30)[PLINK_KEY_COLUMN].nunique() > 1:
        # ".head()" is an optimization so that we don't have to wait for .unique() / .nunique() to run on the full df
        pass
    else:
        if (plink_variants_df[PLINK_KEY_COLUMN].nunique() == 1) and (plink_variants_df[PLINK_KEY_COLUMN].unique()[0] == "."):
            print_non_annotated_plink_file_error()
            exit(30)
    # Figuring out the key
    rsid_column = ask(
        "Which column in PRS file contains the [bold red]rsID[/bold red]?",
        choices=list(pgscatalog_df.columns),
        _config="rsID column",
        _tool=Prompt,
    )
    # Creating keys
    pgscatalog_df[PGSCATALOG_KEY_COLUMN] = pgscatalog_df[rsid_column].astype(str)
    console.print(layout)


def prep_keys_chrompos(pgscatalog_df: pd.DataFrame, plink_variants_df: pd.DataFrame):
    # Figuring out the key
    chrom_column = ask(
        "Which column in PRS file contains the [bold red]chromosome[/bold red]?",
        choices=list(pgscatalog_df.columns),
        _config="chromosome column",
        _tool=Prompt,
    )
    pos_column = ask(
        "Which column in PRS file contains the [bold red]position[/bold red]?",
        choices=list(pgscatalog_df.columns),
        _config="position column",
        _tool=Prompt,
    )
    # Creating keys
    pgscatalog_df[PGSCATALOG_KEY_COLUMN] = pgscatalog_df[chrom_column].astype(str) + ":" + pgscatalog_df[pos_column].astype(str)
    console.print(layout)


def check_the_files_match(pgscatalog_df, plink_variants_df):
    # Heuristics of mapping
    amount_of_shared_keys = len(set(plink_variants_df[PLINK_KEY_COLUMN]).intersection(pgscatalog_df[PGSCATALOG_KEY_COLUMN]))
    amount_of_plink_variants = plink_variants_df.shape[0]
    amount_of_pgscatalog_variants = pgscatalog_df.shape[0]
    # TODO: maybe (!) test whether the *strongest* SNPs are present in the genetic data
    # TODO: to ensure that PRS won't lose it's power?
    if amount_of_shared_keys > 0.50 * amount_of_pgscatalog_variants:
        printout(f"{amount_of_shared_keys} ({amount_of_shared_keys / amount_of_pgscatalog_variants:.0%}) "
                 f"variants from PGSCatalog are available in your plink data ({amount_of_plink_variants} variants)")
        return True
    return False


def maybe_unpack(prs_wm_text_file: str):
    if prs_wm_text_file.endswith(".txt.gz"):
        printout("gzipped PGSCatalog file found, unpacking...")
        process = subprocess.run([
            "gunzip", prs_wm_text_file
        ])
        assert process.returncode == 0, "  ERROR: gunzip crashed"
        prs_wm_text_file = prs_wm_text_file[:-3]
    return prs_wm_text_file


def load_prs_wm(prs_wm_text_file: str):
    # Load, remove # part of the header
    df: pd.DataFrame = pd.read_table(
        prs_wm_text_file,
        sep="\t", engine="c",
        header=0,
        comment='#'  # to skip the PGSCatalog-specific header
    )

    printout(f"  {df.shape[0]} SNPs discovered in the PRS weight matrix")
    if df.shape[0] < 1000:
        printout("  WARNING: This looks like a more classical GRS! But we will continue anyways")

    return df


def annotate_pgscatalog_file(pgscatalog_df):
    # Annotate the essential beta and effect allele columns:
    beta_column = ask(
        "Which column in PRS file contains the [bold red]beta[/bold red]?",
        choices=list(pgscatalog_df.columns),
        _config="beta column",
        _tool=Prompt,
    )
    console.print(layout)

    effect_allele_column = ask(
        "Which column in PRS file contains the [bold red]effect allele[/bold red]?",
        choices=list(pgscatalog_df.columns),
        _config="effect allele column",
        _tool=Prompt,
    )
    console.print(layout)

    columns = {
        beta_column: "effect_weight",
        effect_allele_column: "effect_allele",
    }
    return columns


def render_file_table(df, title):
    df_visual_table = Table(*df.columns, title=title)
    for row in df.values[:5]:
        df_visual_table.add_row(*list(map(str, row)))
    return df_visual_table


def load_data(plink_prefix: str, prs_wm_text_file: str):
    # Unpack, if required
    prs_wm_text_file = maybe_unpack(prs_wm_text_file)

    # Variables
    processed_wm_text_file = prs_wm_text_file + ".annotated"

    # Load PGSCatalog data
    printout("Loading PGSCatalog data and removing comments and indels...")
    pgscatalog_df = load_prs_wm(prs_wm_text_file)
    layout["top"]["leftfile"].update(render_file_table(pgscatalog_df, title="PRS WM file"))
    console.print(layout)
    # Annotate PGSCatalog data
    columns = annotate_pgscatalog_file(pgscatalog_df)
    pgscatalog_df.rename(columns=columns, inplace=True)
    layout["top"]["leftfile"].update(render_file_table(pgscatalog_df, title="PRS WM file"))
    console.print(layout)

    # Load plink data
    printout("Loading plink data...")
    # TODO: load async while waiting for the input above
    plink_variants_df = pd.read_table(
        plink_prefix+".bim",
        sep="\t", header=None,
        # fixed structure of plink19 .bim files - https://www.cog-genomics.org/plink/1.9/formats#bim
        names=["chr_name", "rsid", "dummy_cm_position", "chr_position", "base_allele", "alternative_allele"],
        usecols=["chr_name", PLINK_KEY_COLUMN, "chr_position", "base_allele", "alternative_allele"],
    ) if PLINK_FILEFORMAT[1] == "bim" else pd.read_table(
        plink_prefix+".pvar",
        sep="\t", header=None,
        # fixed structure of plink2 .pvar files - https://www.cog-genomics.org/plink/2.0/formats#pvar
        names=["chr_name", "chr_position", "rsid", "base_allele", "alternative_allele"],
        usecols=["chr_name", "chr_position", PLINK_KEY_COLUMN, "base_allele", "alternative_allele"],
    )

    layout["top"]["rightfile"].update(render_file_table(plink_variants_df, title="plink data"))
    console.print(layout)

    printout("All data are loaded\n-----------------")

    return pgscatalog_df, plink_variants_df, processed_wm_text_file


def preprocess_data(pgscatalog_df, plink_variants_df, processed_wm_text_file):
    # First, ask how to match
    RSID = "rsid"
    CHROMPOS = "pos"
    how_to_match = ask(
        "Would you prefer matching plink data on existing identifiers (rsid) "
        "or create a novel key from chromosome:position pair (pos)?",
        choices=[RSID, CHROMPOS],
        _config="matching on",
        _tool=Prompt,
    )
    if how_to_match == RSID:
        printout("Attempting to match PGSCatalog data with plink data using rsID...")
        prep_keys_rsid(pgscatalog_df, plink_variants_df)
    elif how_to_match == CHROMPOS:
        printout("Attempting to match PGSCatalog data with plink data using chromosome and position...")
        prep_keys_chrompos(pgscatalog_df, plink_variants_df)
    else:
        raise RuntimeError("Unreacheble code")
    was_match_successful = check_the_files_match(pgscatalog_df, plink_variants_df)
    # QC
    # input_files_qc(pgscatalog_df, plink_variants_df)
    # Assess the results and create final matching that plink will be able to use
    pgscatalog_df[[PGSCATALOG_KEY_COLUMN, "effect_allele", "effect_weight"]].to_csv(
        processed_wm_text_file, sep="\t", float_format="%.4e", index=False,
    )
    if not was_match_successful:
        printout("ERROR: Less that a half of the PGSCatalog SNPs are available.")
        printout("ERROR: This may be for one of the following reasons:")
        print_error_files()
        exit(40)


def final_check(pgscatalog_df, plink_variants_df):
    overlaps = pgscatalog_df[PGSCATALOG_KEY_COLUMN].isin(plink_variants_df[PLINK_KEY_COLUMN])
    sample = pgscatalog_df.loc[overlaps, [PGSCATALOG_KEY_COLUMN]].head(10)
    layout["top"]["leftfile"].update(
        Align(
            render_file_table(sample, title="PRS WM file key"),
            align="center",
        )
    )
    layout["top"]["rightfile"].update(
        Align(
            render_file_table(sample, title="plink data key"),
            align="center",
        )
    )
    console.print(layout)
    printout("-----------------\nPlease review the keys the data will be matched on")
    if CONFIG is not None:
        is_to_run = True  # auto mode
    else:
        # this is intentionally not done using "ask()" due to a different logic
        is_to_run = Confirm.ask("Do you confirm you want to match the files on these keys?")  # manual mode
    if is_to_run:
        layout["top"]["leftfile"].update(render_file_table(pgscatalog_df, title="PRS WM file"))
        layout["top"]["rightfile"].update(render_file_table(plink_variants_df, title="plink data"))
        console.print(layout)
        printout("-----------------")
    else:
        printout("You have decided to halt the execution due to the mismatch between the IDs in two files.")
        printout("Please, review this troubleshooting guide, if you need an inspiration about how fixing this:")
        print_error_files()
        exit(50)


def plink_qc(df: pd.DataFrame):
    if PLINK_VERSION >= 2:
        sum_column = "NAMED_ALLELE_DOSAGE_SUM"
    else:
        sum_column = "CNT2"
    avg_snps_used = df[sum_column].mean()
    bad_samples = df.index[df[sum_column] < 0.95 * avg_snps_used]
    if not bad_samples.empty:
        printout("WARNING: The following samples have a lot of missing SNPs "
                 "(more than 5% of SNPs, available on average in other samples)")
        printout("They will dropped from the following analysis")
        printout(",".join(bad_samples))
        df.drop(bad_samples, inplace=True)
    return df


def calculate_prs(plink_prefix: str, processed_wm_text_file: str, output_prefix: str):
    # Variables
    plink_output_file = output_prefix + ".sscore" if PLINK_VERSION >= 2 else output_prefix + ".profile"
    prs_output_file = output_prefix + ".prs"    

    # Running plink
    printout("Running plink PRS scoring...")
    format_arg = "--bfile" if PLINK_FILEFORMAT[0] == "bed" else "--pfile"
    process = subprocess.run([
        PLINK_EXEC, format_arg, plink_prefix, "--score", processed_wm_text_file,
        "1",  # variant ID column index, 1-based
        "2",  # effect allele column index, 1-based
        "3",  # beta column index, 1-based
        "header",
        "--out", output_prefix,
    ], universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if process.returncode != 0 or not os.path.isfile(plink_output_file):
        printout("ERROR: plink has failed")
        printout(process.stderr)
        exit(20)

    # Parsing plink results
    printout("Parsing plink results...")
    try: df = pd.read_table(plink_output_file, sep="\s+", index_col="IID")
    except ValueError: df = pd.read_table(plink_output_file, sep="\s+", index_col="#IID")
    df: pd.DataFrame = plink_qc(df)
    count_column = "ALLELE_CT" if PLINK_VERSION >= 2 else "CNT"
    input_column = "SCORE1_AVG" if PLINK_VERSION >= 2 else "SCORE"
    df[input_column] = df[input_column] * df[count_column]
    scorename = os.path.basename(processed_wm_text_file).split(".")[0]
    output_column = f"SCORESUM_{scorename}"
    df.rename(columns={input_column: output_column}, inplace=True)
    df[[output_column]].to_csv(prs_output_file, sep="\t", float_format="%.4e")
    # Done
    printout(f"PRS is saved to {prs_output_file}")


def printout(text):
    # Logging
    if LOGFILE is not None:
        LOGFILE.write(text+"\n")
    # Printing to interactive console
    # noinspection PyUnresolvedReferences
    old_text = layout["main"].renderable.renderable
    old_lines = old_text.split("\n")
    new_lines = text.split("\n")
    new_text = "\n".join((old_lines + new_lines)[-30:])
    layout["main"].update(Panel(new_text))
    console.print(layout)


def init_render_area():
    # Divide the "screen" in to three parts
    layout.split(
        Layout(name="header", size=3),
        Layout(name="top", size=12),
        Layout(name="main", ratio=1),
        # Layout(name="cmd", size=3),
    )
    # Divide the "side" layout in to two
    layout["top"].split_column(
        Layout(name="leftfile", ratio=1),
        Layout(name="rightfile", ratio=1),
    #    direction="horizontal",
    )

    layout["header"].update(Panel("PRS Tool"))
    layout["main"].update(Panel("Loading..."))
    layout["top"]["leftfile"].update(Panel("PRS WM file will appear here"))
    layout["top"]["rightfile"].update(Panel("Genetic data will appear here"))
    console.print(layout)


def main(args):
    check_input(args)
    plink_prefix, prs_wm_text_file, output_prefix = args.genetic, args.prs_wm, args.out
    # TODO: logging won't work above this line
    global LOGFILE
    LOGFILE = open(output_prefix + ".logging", 'w')
    init_render_area()
    maybe_load_config(args.config)
    check_plink_exec()
    pgscatalog_df, plink_variants_df, processed_wm_text_file = load_data(plink_prefix, prs_wm_text_file)
    # TODO: reuse
    # is_to_reuse = False
    # if os.path.isfile(processed_wm_text_file):
    #     is_to_reuse = Confirm.ask(
    #         "PGS files were preprocessed - would you like to re-use them (recommended) instead of processing again?",
    #     )
    preprocess_data(pgscatalog_df, plink_variants_df, processed_wm_text_file)
    final_check(pgscatalog_df, plink_variants_df)
    calculate_prs(plink_prefix, processed_wm_text_file, output_prefix)
    LOGFILE.close()


if __name__ == "__main__":
    main(arguments)
