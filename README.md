# oskar tools

python tools for working with hdf5 structured datasets obtained using oskar.

## Installation

Tested on Windows 7 using Anaconda and python 2.7.  Should also work with python 3.5.

### prerequisites

positronium (python tools pertaining to positronium)

```
pip install positronium
```

and sspals (python tools for analysing single-shot positron annihilation lifetime spectra)

```
pip install sspals
```

### github.com
```
git clone https://github.com/PositroniumSpectroscopy/oskar_tools
```

and run

```bash
cd ./oskar_tools
python setup.py install
```

This will install the oskar package to your python environment, which can be loaded as usual, e.g.

```python
>>> import oskar
```

setup.py also installs console scripts for configuring oskar and analysing data for a given run.  These are

|script         | description                              |
|---------------|------------------------------------------|
|oskar_dset     | configure defaults.json                  |
|oskar_info     | print info about a run                   |
|oskar_average  | average DataFrame data                   |
|oskar_vrange   | find the vertical range of waveform data |
|oskar_sspals   | sspals analysis of waveform data         |
|oskar_count    | count trigger events in waveform data    |

These should be accessible from a command prompt.  For example

```bash
oskar_info -r 20160417_164136

  20160417_164136
     Author:       AA
     Description:  Image positron beam on MCP
  SQUIDS: 4
  VARS: [u'HOLD']
  REC: []
  Unique VAR combinations: 2
  Number of loops: 2.0
```

Run with option --help for more information on usage.

To run these scripts inside of a Jupyter notebook see 'examples/Run Analysis'

---

## defaults.json

Default configurations for running scripts can be stored in the defaults.json file, which is stored in the oskar_tools directory.  The most important option is 'base' as this specifies the base directory in which data is stored.

default.json can be created manually, however, it's simplest to modify it is using the command line tool oskar_dset, which should be installed to your system path when running setup.py.  For instance, to write the 
attribute "base" with value "Z:\\\\" and to print the result, open a cmd prompt and type

```
oskar_dset -a "base" "Z:\\Data" --pprint
```

Also, see the notebook 'examples/Defaults class'.

Example json file:

```json
{
     "count": [ 
          "CH_A0",
          "CH_A1"
     ],
     "average": [
          "{'SSPALS':['t0','DF','Range','FWHM']}"
     ],
     "base": [
          "Z:\\Data"
     ],
     "rid": [
          "20160224_192014"
     ]
}
```