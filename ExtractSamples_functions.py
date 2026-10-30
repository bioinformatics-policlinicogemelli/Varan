#Copyright 2025 bioinformatics-policlinicogemelli

#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

"""Extract (-e): sample list checks and the copy of a study's data files
restricted to the requested samples. The per-file rules live in
study_files.py, shared with remove and update.
"""

import sys
from pathlib import Path

import pandas as pd
from loguru import logger

from study_files import subset_study



def check_sample_list(extract_path: str, oldpath: str) -> None:
    """Validate the sample list file and check consistency with clinical sample data.

    This function performs several checks to ensure the list of samples to be extracted.
    It checks if this list contains only a single column; Includes at least one sample
    present in the 'data_clinical_sample.txt' file; Warns if some of the samples in the
    list are not found in the clinical sample data; Warns if the list includes all
    samples from the original study.

    Args:
        extract_path (str): Path to the text file containing sample ID to extract.
                            The file should contain one sample ID per line.
        oldpath (str): Path to the directory containing the original study data,
                            including the 'data_clinical_sample.txt' file.

    Returns:
        None

    """
    with Path(extract_path).open() as sample_list:
        first_line = sample_list.readline()
        if len(first_line.split("\t")) > 1:
            logger.critical(
                f"The file {extract_path} contains more than a column."
                f"It may not be in the correct format!")
            sys.exit(1)

    with Path(extract_path).open() as sample_list:
        all_samples_to_extract = {sample.strip() for sample in sample_list}

        sample = pd.read_csv(
            (Path(oldpath) / "data_clinical_sample.txt"), sep="\t", skiprows=4)
        old_samples = set(sample["SAMPLE_ID"])

        if not any(sample in old_samples for sample in all_samples_to_extract):
            logger.critical(
                "The sample(s) you are trying to extract are not"
                "present in data_clinical_sample.txt file! Please check again!")
            sys.exit(1)

        missing_samples = set(all_samples_to_extract) - set(old_samples)
        if missing_samples:
            logger.warning(
                "Some of the samples you are trying to extract"
                "may not be present in data_clinical_sample.txt file!")
            logger.warning(f"Missing sample(s): {', '.join(missing_samples)}")

        if len(set(old_samples) - set(all_samples_to_extract)) == 0:
            logger.warning(
                "It looks like you are extracting all the samples from"
                "the original study, but this might be redundant, as it"
                "replicates the existing one.")


def extract_all_data(oldpath: str, sample_ids: list, output: str) -> None:
    """Write every per-sample data file of the study at `oldpath` into
    `output`, keeping only `sample_ids` - see study_files.subset_study()."""
    subset_study(oldpath, sample_ids, output, keep=True)
