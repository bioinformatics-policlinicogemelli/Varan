#cbioportal validator
git clone https://github.com/cBioPortal/cbioportal-core.git
mv cbioportal-core/scripts/importer .
rm -r cbioportal-core

#oncokb annotator
git clone https://github.com/oncokb/oncokb-annotator.git
pip install -r oncokb-annotator/requirements/pip3.txt
pip install -r oncokb-annotator/requirements/common.txt

#requirements for Varan
pip install -r requirements.txt