# HII Region Database Generator
Generates a database containing all relevant information about Galactic HII regions in the
[WISE Catalog of Galactic HII Regions](http://astro.phys.wvu.edu/wise/)
([paper](https://ui.adsabs.harvard.edu/abs/2014ApJS..212....1A/abstract)) and recent
radio recombination line surveys. The database can be [downloaded through the
Harvard Dataverse](https://doi.org/10.7910/DVN/NQVFLE). The structure of the database
is described in `schema.txt`.

## Requirements & Installation
The following packages are required for this package to work "out-of-the-box"
1. `astropy`
2. `numpy`
3. `kd` (from https://github.com/tvwenger/kd)

The easiest way to run this code is
```bash
pip install git+https://github.com/tvwenger/pyqt-fit.git
pip install git+https://github.com/tvwenger/kd.git
git clone https://github.com/tvwenger/hii_db.git
```

## Usage
This repository includes a CSV version of the WISE Catalog of Galactic HII Regions, which
can be used to populate a "WISE-only" version of the database (i.e., a database that does
not include any ancillary radio recombination line data). This is accomplished via
`generate.py` with the `--wiseonly` flag.

```bash
python generate.py --help
```
```
usage: generate.py [-h] [--wise WISE] [--wiseonly] db

HII Region Database Generator

positional arguments:
  db           The database filename

optional arguments:
  -h, --help   show this help message and exit
  --wise WISE  WISE Catalog CSV filename (default: wise/wise_hii_V3.0_hrds.csv)
  --wiseonly   Generate database with WISE Catalog only (default: False)
```
```bash
python generate.py new_database.db --wise wise/wise_hii_V3.0_hrds.csv --wiseonly
```

To add the parallax data from Reid et al. (2019) to the (existing) database, use `add_parallax.py`:
```bash
python add_parallax.py --help
```
```
usage: add_parallax.py [-h] [--data DATA] [--refs REFS] db

HII Region Database Parallax Table Generator

positional arguments:
  db           The database filename

optional arguments:
  -h, --help   show this help message and exit
  --data DATA  The parallax data filename (default: data/reid_2019/reid2019_merge.txt)
  --refs REFS  The parallax data references filename (default: data/reid_2019/reid2019_refs.txt)
```
```bash
python add_parallax.py new_database.db --data data/reid_2019/reid2019_merge.txt --refs data/reid_2019/reid2019_refs.txt
```

To add Wenger et al. (2018) Monte Carlo kinematic distances to the (existing) database, use `add_distances.py`:
```bash
python add_distances.py --help
```
```
usage: add_distances.py [-h] [-n NUM_SAMPLES] [-r ROTCURVE] [-t TABLENAME] db

HII Region Database Distances Table Generator

positional arguments:
  db                    Database filename

optional arguments:
  -h, --help            show this help message and exit
  -n NUM_SAMPLES, --num_samples NUM_SAMPLES
                        Number of Monte Carlo samples (default: 5000)
  -r ROTCURVE, --rotcurve ROTCURVE
                        Rotation curve (default: reid19_rotcurve)
  -t TABLENAME, --tablename TABLENAME
                        Table name (default: Distances_Reid2019)
```
```bash
python add_distances.py new_database.db -n 5000 -r reid19_rotcurve -t Distances_Reid2019
```

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
