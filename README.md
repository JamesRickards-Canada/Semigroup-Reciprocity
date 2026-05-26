# Semigroup-Reciprocity

This repository is supplementary to the paper [*Reciprocity obstructions in semigroup orbits in SL(2, Z)*](https://doi.org/10.1215/00127094-2025-0017) ([Arxiv](https://arxiv.org/abs/2401.01860)), by James Rickards and Katherine E. Stange.

The goal of the package is to efficiently compute data about certain orbits of semigroups, including computing:
* All matrices in the semigroup up to a certain bound;
* The orbit of a semigroup on a vector up to a certain bound;
* Estimation of the growth rate of the semigroup orbits;
* The set of missing integers in a semigroup orbit.

See [here](USERGUIDE.md) for a more detailed Users Guide.

## Bibtex

If you use this code in a project, please let me know! Here is a Bibtex entry that points to the latest version:

## Installation Instructions

### Downloading the code
Call ```git clone https://github.com/JamesRickards-Canada/Semigroup-Reciprocity.git```. If you are on Windows, be sure to git clone from WSL (see below), as Windows line endings (carriage returns) may be added to files, causing issues.

### Prerequisites
* PARI/GP, but _not_ the downloaded ready-to-go binary. The PARI/GP website has binaries for Windows and Mac avaliable, but these will not work with the package. See below for OS specific instructions.
* You should have a guess as to the location of the ```pari.cfg``` file for the version of PARI/GP you are running. Suggestions on how to do this can be found below.
* This package has been tested on versions:
<p align="center"> 2.15.3 - 2.15.5, &emsp; 2.17.0 - 2.17.3, &emsp; 2.18.1</p>
If you try to use it on a different version, be aware that it may not work.

### Operating systems
* **Linux** - No further requirements
* **Windows** - You need to use Windows Subsytem for Linux. See the [guide](https://pari.math.u-bordeaux.fr/PDF/PARIGP-Windows.pdf) I wrote for additional instructions.
* **Mac** - You need to have [Homebrew](https://brew.sh/) installed. This is also an easy way to install PARI/GP: ```brew install pari```

### Configuring and building the package
* From inside the project folder, call ```./configure``` to initialize the project. This helps you search for ```pari.cfg```, and stores the location to a file. You should supply it with a folder to search in!
* The script displays the corresponding versions of the found files, so if you have multiple versions, you can choose the correct one. This can be useful if you keep multiple copies of PARI/GP around.
* If the location of the installation of PARI/GP does not change, you do not need to reconfigure. If when you update PARI/GP there is a new location (e.g. if the version number is in the file path of ```pari.cfg```), you should call ```./configure``` again.
* Call ```make``` to build the project, and ```make clean``` to remove all .o object files. If you change versions of PARI/GP, you must ```make clean``` and then ```make``` again.
* Once this is done, a call to ```gp semigroup``` starts gp with the package installed!
* Further updates to the package (that you download with ```git pull```) just require calling ```make``` again.
* Call ```?semigroup``` to access further help for this package.

### Where is pari.cfg?
* The configuration file will search for this, but it is preferrable to not search your entire hard drive (as this can be very slow). So, you should at least supply a guess as to the location of ```pari.cfg```. Often only the top-level folder (e.g. ```/usr``` or ```/opt```) suffices.
* On Linux or WSL, if you build PARI/GP from source, it should be located in ```/usr/local/lib/pari/pari.cfg```, or at least somewhere in the ```/usr``` folder.
* On a Mac, if you install PARI/GP with Homebrew, it may be found in a folder like ```/opt/homebrew/Cellar/pari/VERSION/lib/pari```. Searching ```/opt/homebrew``` should be fine.
* If you are obtaining it through SageMath, it might be found where the library files of SageMath are
* Assuming you open PARI/GP with the command ```gp```, try ```type -a gp```, which will display where this command lives. The corresponding file(s) are likely symbolic links, and you can call ```readlink -f LOCATION``` on each of them to see where it lives. This can provide a clue as to the place to search for ```pari.cfg```.
* Another clue comes from gp itself. Open gp, and type ```default()```. Look for the entries ```datadir``` and ```help```. It is _often_ the case that ```datadir``` is in ```X/share/pari```, ```help``` is in ```X/bin/gphelp```, and ```pari.cfg``` is in ```X/lib/pari/pari.cfg```.

### Troubleshooting

* When in doubt, ```make clean``` to reset the object files, reconfigure with ```./configure```, and then run ```make``` again.
* If you have any further issues with installation or using the code, please let me know!
 
