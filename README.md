# At-the-time and Back-in-time Persistent Sketches

This anonymous repo includes the experiment code for ATTP and BITP sketches.

## Dependencies
c/c++ compileres: must be gcc/g++ >= 8. Sorry, the code won't work with other compilers.
libraries: lapacke, cblas, fftw3

## How to compile:
1. To configure it with -O2 -DNDEBUG flags, run

    ./configure
    
or to configure it with -O0 -g flags, run

    ./configure --enable-debug 

2. To compile, run

    make

3. If you added/removed some files and/or added headers to some cpp, call the following first before compiling it

    make depend

## How to obtain datasets

For heavy hitters, we use the 1998 world-cup dataset. For matrix sketches, we
use synthetic datasets.

### 1998 World-cup server logs

Available [here](ftp://ita.ee.lbl.gov/html/contrib/WorldCup.html).

See data\_proc/world-cup/prepare\_data.sh for how we preprocessed the logs.

Generated files:

    - data/world-cup-98-client-id.txt: (TIMESTAMP, Client-ID) pairs

    - data/world-cup-98-object-id.txt: (TIMESTAMP, Object-ID) pairs

    - data/wc-client-id-attp-w-stats-report.txt: Client-ID data and queries
      for ATTP.

    - data/wc-client-id-bitp-w-stats-report.txt: Client-ID data and queries
      for BITP.

    - data/wc-object-id-attp-w-stats-report.txt: Object-ID data and queries
      for ATTP.

    - data/wc-object-id-bitp-w-stats-report.txt: Object-ID data and queries
      for BITP.

### matrix sketch dataset

run data\_proc/gen\_mat\_data.sh

Generated files:
    
    - data/matrix\_small-w-stats-report: d = 100, n = 50000

    - data/matrix\_medium-w-stats-report: d = 1000, n = 50000

    - data/matrix\_big-w-stats-report: d = 10000, n = 50000, don't run exact query on this dataset or you'll need a very large memory

## How to run a query

    ./driver run configs/xxx.conf

To get help, run

    ./driver

or
    
    ./driver help <QueryType>

