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

"""Remove (-r): sample list checks and the copy of a study's data files
without the given samples. The per-file rules live in study_files.py,
shared with extract and update.
"""

import sys
from pathlib import Path

import pandas as pd
from loguru import logger

from study_files import subset_study



def check_sample_list(remove_path: str, oldpath: str) -> None:
    """Validate that sample list is correctly formatted and exists in the dataset.

    Args:
        remove_path (str): Path to file listing sample IDs to remove.
        oldpath (str): Directory containing the clinical sample file.

    Returns:
        None

    Raises:
        SystemExit: If the sample list format is invalid or results in an empty study.

    """
    with Path(remove_path).open(encoding="utf-8") as sample_list:
        first_line = sample_list.readline()
        if len(first_line.split("\t")) > 1:
            logger.critical(
                f"The file {remove_path} contains more than a column. "
                "It may not be in the correct format!")
            sys.exit(1)

    with Path(remove_path).open(encoding="utf-8") as sample_list:
        all_samples_to_remove = {sample.strip() for sample in sample_list}

        sample_path = Path(oldpath) / "data_clinical_sample.txt"
        sample = pd.read_csv(sample_path, sep="\t", skiprows=4)
        old_samples = set(sample["SAMPLE_ID"])

        if not any(sample in old_samples for sample in all_samples_to_remove):
            logger.critical("The sample(s) you are trying to remove are not present "
            "in data_clinical_sample.txt file! Please check again!")
            sys.exit(1)

        missing_samples = set(all_samples_to_remove) - set(old_samples)
        if missing_samples:
            logger.warning("Some of the samples you are trying to remove may not be "
            "present in data_clinical_sample.txt file!")
            logger.warning(f"Missing sample(s): {', '.join(missing_samples)}")

        if len(set(old_samples) - set(all_samples_to_remove)) == 0:
            logger.critical("It looks like you are removing all the samples from "
            "original study! Cannot create an empty study!")
            sys.exit(1)


def delete_all_data(oldpath: str, sample_ids: list, output: str) -> None:
    """Write every per-sample data file of the study at `oldpath` into
    `output`, without `sample_ids` - see study_files.subset_study()."""
    subset_study(oldpath, sample_ids, output, keep=False)
