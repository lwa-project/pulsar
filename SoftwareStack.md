# The LWA Pulsar Stack

Pulsar data reduction software (and versions) currently avaliable on the LWA Users Computing Facility ([LWA Memo #193](http://www.phys.unm.edu/~lwa/memos/memo/lwa0193d.pdf)).  This is in addition to the [LWA Software Library](https://fornax.phys.unm.edu/lwa/trac/) and the [pulsar extension](https://github.com/lwa-project/pulsar/).

## TEMPO
http://tempo.sourceforge.net/
```
git clone git://git.code.sf.net/p/tempo/tempo
git checkout 61c8011b0477f005dd4f05ad037fae83c51eafb4
./prepare
./configure
make
sudo make install
```

## PRESTO
https://www.cv.nrao.edu/~sransom/presto/
```
git clone https://github.com/scottransom/presto.git
git checkout v4.0
cd src
make
```

## psrfits_utils
```
git clone https://github.com/lwa-project/psrfits_utils.git
git checkout d26cac580477a2f2148e1a69fe70178747bf7ed1
./prepare
PYTHON_VERSION="3.8" ./configure
make
sudo make install
```

## EPSIC
```
git clone https://github.com/straten/epsic.git
git checkout 5315cc634f6539ea0a34e403e492472b97e0f086
cd src/
./bootstrap
./configure
make
sudo make install
```

## PSRCAT
https://www.atnf.csiro.au/research/pulsar/psrcat/
```
wget https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz
tar xzvf psrcat_pkg.tar.gz
mv psrcat_tar psrcat
./makeit
sudo cp psrcat /usr/local/bin/
```

## PSRCHIVE
http://psrchive.sourceforge.net/
```
git clone git://git.code.sf.net/p/psrchive/code
git checkout a1709aab0948935dcaf5305c0aefd129cf954f8b
unset TEMPO2
./bootstrap
PYTHON=/usr/bin/python3.8 ./configure --enable-shared --disable-tempo2 --with-psrcat=/usr/local/psrcat/
make
sudo make install
```

## DSPSR
http://dspsr.sourceforge.net/
```
git clone git://git.code.sf.net/p/dspsr/code
git checkout 6c762e6fc6bc2a1595000f34bcfc0ca0c43e9954
unset TEMPO2
sed -e 's/mwa/lwa mwa/g' -i ./config/backends.default 
./bootstrap
./configure --enable-shared --with-cuda-dir=/usr/local/cuda --with-cuda-include-dir=/usr/local/cuda/include --with-cuda-lib-dir=/usr/local/cuda/lib64
make
sudo make install
```
