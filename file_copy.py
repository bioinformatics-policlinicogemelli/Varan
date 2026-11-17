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

"""Copy files into an output study without carrying over their permissions.

shutil.copy and shutil.copy2 also copy the permission bits, so a file
taken from a read-only study (or from a read-only installation of Varan)
stays read-only in the output, and Varan can no longer rewrite it there.
copy_content copies only the data: the new file gets the default
permissions of the user running Varan.
"""

import shutil
from pathlib import Path


def copy_content(src, dst) -> Path:
    """Copy the content of src to dst (a file path or an existing folder).

    Returns the path of the copied file.
    """
    dst = Path(dst)
    if dst.is_dir():
        dst = dst / Path(src).name
    shutil.copyfile(src, dst)
    return dst
