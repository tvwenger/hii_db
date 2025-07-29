# HII Region Database Generator
Generates a database containing all relevant information about Galactic HII regions in the
[WISE Catalog of Galactic HII Regions](http://astro.phys.wvu.edu/wise/)
([paper](https://ui.adsabs.harvard.edu/abs/2014ApJS..212....1A/abstract)) and recent
radio recombination line surveys. The database can be [downloaded through the
Harvard Dataverse](https://doi.org/10.7910/DVN/NQVFLE). The structure of the database
is described in `schema.txt`.

## Data Source
The data needed to generate the database can be downloaded from the [Harvard Dataverse]().

## Installation
Run this package in a new `conda` environment:
```bash
mamba env create -n hii_db -c conda-forge "python<3.12" "numpy<2.0.0" scipy matplotlib pip pandas jupyter astropy
mamba activate hii_db
pip install git+https://github.com/tvwenger/pyqt-fit.git#egg=pyqt-fit[cython]
pip install git+https://github.com/tvwenger/kd.git
```

Then clone this repository:
```bash
git clone https://github.com/tvwenger/hii_db.git
```

## Usage
The database is built interactively with the notebook `generate_database.ipynb`. See that notebook for more information.

## Issues and Contributing

Please submit issues or contribute to the development via [Github](https://github.com/tvwenger/hii_db).

## License and Warranty

GNU Public License
http://www.gnu.org/licenses/

This package is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This package is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this package. If not, see http://www.gnu.org/licenses/
