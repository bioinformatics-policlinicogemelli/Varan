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

"""Varan desktop app: a window that runs Varan's Docker image.

Runs on the user's computer (Windows, macOS, Linux) with nothing but the
Python standard library and tkinter. Varan itself, VEP and vcf2maf stay in
the Docker image the command line uses, so both give the same results:
this package only collects the inputs, translates the paths for the
container, starts it and shows its log and report.

    python -m varan_gui
"""

__version__ = "0.1.0"
