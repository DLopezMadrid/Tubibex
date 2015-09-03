# TUBIBEX v1.2
### Developed by Daniel Lopez Madrid

Tubibex is an algorithm for calculating capture tubes based on the ibex library

My MSc in Advanced Control and Systems engineering is in the doc folder:
D. Lopez, Development and implementation of an algorithm for the calculation of positive invariant sets using interval methods, University of Manchester, 2015

Capture tube paper: L. Jaulin, D. Lopez, Le Doze, S. Le Menec, J. Ninin, G. Chabert, M. S. Ibnseddik, A. Stancu (2015), Computing capture tubes, Reliable Computing, accepted. PDF: https://www.ensta-bretagne.fr/jaulin/paper_captube.pdf



Available in /doc folder


Ibex library: 	http://www.ibex-lib.org


## INSTALL:
Tubibex was developed in Linux Mint 17.2 x64 and it uses the following software packages:
	
*  QT 4.8.6
    * https://download.qt.io/archive/qt/4.8/4.8.6/
	* http://doc.qt.io/qt-4.8/index.html

* gdb 7.9.1
    * https://www.sourceware.org/gdb/download/
		
    * install commands (needs to be python-enabled)
      * ./configure --prefix=/usr --with-python
      * make
      * make install

* Ibex 2.1.16
  * http://www.ibex-lib.org/doc/install.html#standard-install-stable

  * The following applications must be installed.
    * g++
    * gcc
    * flex
    * bison
    * python2.x (warning: the script are currently not compatible with python3)
    * make
    * pkg-config (optionnal)
    * jdk (optionnal)
    * On Ubuntu, you can install all you need with:
      * ~$ sudo apt-get install -y python2.7 flex bison gcc g++ make pkg-config


* VIBes 0.2.0 beta
  * https://github.com/ENSTABretagneRobotics/VIBES/releases
  * to compile the viewer in linux use cmake 
  * cmake ./CMakeLists.txt
  * make

* DynIbex for Ibex 2.1.16
    * http://perso.ensta-paristech.fr/~chapoutot/dynibex/


* Substitute the ibex_IntervalVector.h and ibex_simulation.h files with the ones in tubibex/3rd_party/ibex_mod/
the files should be in the following location:
Windows: C:\MinGW\msys\1.0\home\[USERNAME]\Ibex\ibex-2.1.16\include\ibex
Linux: /usr/local/include/ibex

GPL v2.0 LICENSE
