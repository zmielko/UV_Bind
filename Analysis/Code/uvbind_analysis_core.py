"""Functions for UV Bind analysis.

Functions used across different ipython notebook files in the UV Bind
analysis.



"""

from __future__ import annotations
from collections import namedtuple
import math

from bokeh.io import export_svg
from bokeh.plotting import figure, show, ColumnDataSource
from bokeh.models import CDSView, GroupFilter
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.sandbox.regression.predstd import wls_prediction_std

PRIMER = "CGCACACATACACATACACACATACACATACGCAC"
PYDI = ("TT", "TC", "CT", "CC", "AA", "GA", "AG", "GG")

Pydi_Tuple = namedtuple("Pyrimidine_Dinucleotides", ["TT", "TC", "CT", "CC"])
Pydi_Tuple_rc = namedtuple("Pyrimidine_Dinucleotides", ["AA", "GA", "AG", "GG"])


def reverse_complement(dna_sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complements = {"T": "A", "A": "T", "C": "G", "G": "C"}
    reverse = dna_sequence[::-1]
    result = ''
    for i in reverse:
        result += complements[i]
    return result


def count_dna_substring(string: str, substring: str) -> int:
    """Count occurances of a substring and its reverse complement."""
    count = 0
    substring_rc = reverse_complement(substring)
    for i in range(len(string) - len(substring) + 1):
        if string[i:i+len(substring)] in (substring, substring_rc):
            count += 1
    return count


def count_overlapping_kmers(string: str,
                            kmer_set: set,
                            k: int) -> int:
    """Count k-mer occurances in a string from a set."""
    count = 0
    for i in range(len(string) - k + 1):
        query = string[i:i+k]
        if query in kmer_set:
            count += 1
    return count


def read_kmer_file_pair(file_a: str,
                        file_b: str,
                        suffixes: tuple[str, str]) -> pd.DataFrame:
    """Return a merged dataframe from a tuple of 2 kmer files.

    K-mer files are defined as a tab seperated file with the first column as
    the k-mers in one orientaion and the second as the other. This takes two
    files, performs an inner join on the k-mer columns and returns a merged
    dataframe with defined suffix labels for both files. The column labels
    are based on k, so 8-mer, 7-mer, or 6-mer.
    """
    dataframe_a = pd.read_csv(file_a, sep='\t')
    dataframe_b = pd.read_csv(file_b, sep='\t')
    merge_on = list(dataframe_a.columns[:2])
    result = pd.merge(dataframe_a, dataframe_b, on=merge_on, suffixes=suffixes)
    return result


def plot_range_from_x_y(x: list[float],
                        y: list[float],
                        int_range: bool = True,
                        extend_percent: float = 0.01) -> tuple[float, float]:
    """Return a tuple for a plot draw range.

    Given a dataframe and two column names, returns a tuple with
    a min and max value to draw a scatterplot of the data in
    the two columns. The range is the min and max of all values
    extended by the percentage from extend_percent.
    """
    max_value = max(((max(x), max(y))))
    min_value = min(((min(x), min(y))))
    range_value = max_value - min_value
    extend = range_value * extend_percent
    if int_range:
        plot_range = (math.floor(min_value - extend),
                      math.ceil(max_value + extend))
    else:
        plot_range = (min_value - extend, max_value + extend)
    return plot_range


def filter_alldata_dataframe(dataframe: pd.DataFrame,
                             primer: str = PRIMER) -> pd.DataFrame:
    """Remove flagged sequence and non-double stranded sequences.

    Given a dataframe with column names in the alldata format, returns
    a dataframe with all flagged rows, rows with no sequence, and
    rows with a sequence without a primer removed.
    """
    alexa = dataframe.columns[7]
    # Filter out flags
    dataframe = dataframe[(dataframe["Cy3Flags"].isnull()) &
                          (dataframe[f"{alexa}Flags"].isnull())]
    # Remove rows with no sequence
    dataframe = dataframe[dataframe["Sequence"].notnull()]
    # Remove rows without primer
    dataframe = dataframe[dataframe["Sequence"].str.endswith(primer)]
    dataframe = dataframe.reset_index(drop=True)
    return dataframe


def has_pydi(string: str,
             pydi_tuple=PYDI) -> bool:
    """Return True if string contains a pyrimidine dinucleotide else false."""
    pydi_set = set(pydi_tuple)
    for i in range(len(string) - 1):
        if string[i:i+2] in pydi_set:
            return True
    return False


def replicate_values(dataframe: pd.DataFrame,
                     column: str,
                     min_count: int = 2) -> pd.DataFrame:
    """Return a dataframe of values with a minimum number of replicates.

    Given a dataframe and column, finds the number of replicates of that
    column and returns a dataframe with a minimum amount of replicates.
    """
    # Calculate value counts and return as a dataframe
    count = pd.DataFrame(dataframe[column].value_counts()).reset_index()
    # Rename the columns
    count = count.rename(columns={column: "Count"})
    count = count.rename(columns={"index": column})
    # Query for items that have a minimum value count
    count = count.query("Count >= @min_count")
    return count


def process_alldata_file(file_path: str,
                         aggregate_medians: bool = True,
                         ln_transform: bool = True,
                         label_has_pydi: bool = True) -> pd.DataFrame:
    """Given an alldata file, return a dataframe of custom sequences.

    Takes an alldata file and does the following pipeline:
    1. Read as a dataframe
    2. Remove universal sequences
    3. Remove trailing replicate labels so replicate data can be queried
    4. Subset the column into the sequence name, sequence, and signal values
    5. Calculate median values of the replicates
    6. Transform into natural log space
    7. Add a column indicating if there is a pyrimidine dinucleotide

    Returns a dataframe with a Name, Sequence, Signal, Has_PyDi columns.
    """
    # Read dataframe
    dataframe = pd.read_csv(file_path, sep='\t')
    alexa = dataframe.columns[7]
    dataframe = filter_alldata_dataframe(dataframe)
    # Remove universal UV-Bind Probes
    dataframe = dataframe[~dataframe["Name"].str.startswith("All")]
    dataframe = dataframe.reset_index(drop=True)
    # Remove trailing replicate label for median calculation
    dataframe["Name_Group"] = dataframe["Name"]\
        .apply(lambda x: '_'.join(x.split('_')[:-1]))
    # Filter for sequences with 2+ replicates
    replicate_counts = replicate_values(dataframe, "Name_Group")
    name_set = set(replicate_counts["Name_Group"])
    dataframe = dataframe[dataframe["Name_Group"]
                          .isin(name_set)].reset_index(drop=True)
    if aggregate_medians:
        dataframe["Name"] = dataframe["Name_Group"]
        del dataframe["Name_Group"]
    # Subset by name, sequence and signal columns
    dataframe = dataframe[["Name", "Sequence", alexa]]
    dataframe = dataframe.rename(columns={alexa: "Signal"})
    # Aggregate into medians
    if aggregate_medians:
        dataframe = dataframe.groupby(by=["Name", "Sequence"]).aggregate(np.median)
        dataframe = dataframe.reset_index()
    else:
        dataframe = dataframe.reset_index(drop=True)
    # Transform into natural log space
    if ln_transform:
        dataframe["Signal"] = dataframe["Signal"].apply(lambda x: np.log(x))
    # Add column indicating if a sequence can form a dimer
    if label_has_pydi:
        dataframe["Has_PyDi"] = dataframe["Sequence"].apply(lambda x: has_pydi(x))
    return dataframe


def scale_column(dataframe: pd.DataFrame,
                 column_name: str,
                 slope: float,
                 intercept: float) -> pd.DataFrame:
    """Scale a column in a dataframe.

    Given a dataframe, column_name, slope, and intercept, scales a column's
    values such that the slope is equal to 1.
    """
    # Save the pre-scale values in another column
    dataframe[f"{column_name}_PreScale"] = dataframe[column_name]
    # Apply the transformation
    dataframe[column_name] = dataframe[column_name]\
        .apply(lambda x: (x - intercept) / (slope))
    dataframe = dataframe.reset_index(drop=True)
    return dataframe


def scale_uv_on_non_uv(non_uv_file, uv_file):
    """From a non-UV and UV file, scales the UV values and return the results.

    Given a non-UV and UV file, performs the following pipeline:
    1. Read and process the files as done in process_alldata_file
    2. Calculate the median values of replicates
    3. Merge them into a single dataframe
    4. Perform a linear regression for sequences that cannot form
        pyrimidine dimers.
    5. Scale the UV values so that the slope of the regresion is 1.

    Returns a dataframe with both prescale and scaled UV values.
    """
    # Read data for Non_UV and UV
    non_uv_df = process_alldata_file(non_uv_file)
    uv_df = process_alldata_file(uv_file)
    # Combine dataframes and do regression of non-damageable probes
    merged_df = pd.merge(non_uv_df,
                         uv_df,
                         on=["Name", "Sequence", "Has_PyDi"],
                         suffixes=("_Non_UV", "_UV"))
    # Perform linear regression
    merged_df_no_pydi = merged_df[~merged_df["Has_PyDi"]]
    regression = stats.linregress(merged_df_no_pydi["Signal_Non_UV"],
                                  merged_df_no_pydi["Signal_UV"])
    # Scale dataframe
    merged_df = scale_column(merged_df,
                             "Signal_UV",
                             regression.slope,
                             regression.intercept)
    return merged_df


def ols(dataframe: pd.DataFrame,
        column_x: str,
        column_y: str,
        alpha: float,
        column_damage: str) -> pd.DataFrame:
    """Ordinary least squares regression.

    Given a dataframe, x values, y values, and an alpha, runs an ols model
    and returns a copy of the dataframe with the following new colums:
    1. Prediction Interval Upper Bounds
    2. Prediction Interval Lower Bounds
    3. Predicted values

    The model is trained on columns x and y for rows where those in
    column_damage are equal to 0.
    """
    # Query for non-damageable sequences
    copy_df = dataframe.copy()
    copy_df = copy_df[~copy_df["Has_PyDi"]]
    # Fit OLS model on the data
    ols_model = smf.ols(f"{column_y} ~ {column_x}",
                        data=copy_df).fit(cov="HC3")
    # Calculate a prediction interval
    pi_std, pi_lower, pi_upper = wls_prediction_std(ols_model, alpha=alpha)
    # Extend the prediction to the whole dataset range
    e_input = dataframe.copy()
    e_input["Weight"] = 1
    e_input = e_input[["Weight", column_x]]
    pi_std, pred_lower, pred_upper = wls_prediction_std(res=ols_model,
                                                        exog=e_input,
                                                        alpha=alpha)
    prediction_array = ols_model.predict(exog=e_input)
    # Add additional columns of the interval and predicted values
    dataframe["Prediction_Upper"] = pred_upper
    dataframe["Prediction_Lower"] = pred_lower
    dataframe["Predicted"] = prediction_array
    return dataframe


def hamming_distance(string_a: str, string_b: str) -> int:
    """Calculate the hamming distance between two strings."""
    distance = 0
    for idx, i in enumerate(string_a):
        if string_b[idx] != i:
            distance += 1
    return distance


def unique_kmers(string: str, k: int) -> set[str]:
    """Given a string and k, return all unique kmers in a string as a set."""
    kmer_set = set()
    for i in range(len(string) - k + 1):
        kmer_set.add(string[i:i+k])
    return kmer_set


def scatterplot_of_groups(plot_df,
                          plot_range,
                          circle_size,
                          palette,
                          ticker,
                          groups,
                          draw_pred_int=False,
                          pred_column_x="Signal_Non_UV",
                          pred_upper="Prediction_Upper",
                          pred_lower="Prediction_Lower"):
    """"""
    # Create bokeh plot object
    p = figure(width=800, height=800, x_range=plot_range, y_range=plot_range)
    source = ColumnDataSource(plot_df)
    # For each group and their corresponding palette color, draw datapoints
    for i, color in zip(groups, palette):
        view = CDSView(source=source,
                      filters=[GroupFilter(column_name="Group", group=i)])
        p.circle("Signal_Non_UV",
                 "Signal_UV",
                 source=source,
                 view=view,
                 size=circle_size,
                 color=color)
    # Draw lines for a prediction interval if specified
    if draw_pred_int:
        line_df = plot_df.sort_values(by=pred_column_x).reset_index(drop=True)
        for i in (pred_upper, pred_lower):
               p.line(line_df[pred_column_x],
               line_df[i],
               line_dash="dashed",
               color="black",
               line_width=3)
    p.xaxis.ticker = ticker
    p.yaxis.ticker = ticker
    p.xaxis.major_label_text_font_size = '0pt'
    p.yaxis.major_label_text_font_size = '0pt'
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.xaxis.minor_tick_line_color = None
    p.yaxis.minor_tick_line_color = None
    p.toolbar_location = None
    return p
