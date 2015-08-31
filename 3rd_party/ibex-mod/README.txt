INSTALLATION INSTRUCTIONS FOR WINDOWS & LINUX


1-Install ibex-2.1.16 and soplex-1.7.2 as explained here:
http://www.ibex-lib.org/doc/install.html

2-Substitute the file ibex_IntervalVector.hm and ibex_simulation.h with the one attached.
If you have followed the steps in the install guide above, the file should be in the following location:
Windows: C:\MinGW\msys\1.0\home\[USERNAME]\Ibex\ibex-2.1.16\include\ibex
Linux: /usr/local/include/ibex

3-Open the bubbibex project using QT (tested with versions 4.7.0 and 4.8.6) 
Download: https://download.qt.io/archive/qt/4.8/4.8.6/
Installation: http://doc.qt.io/qt-4.8/installation.html

4-Make sure that the bubbibex.pro file is using the correct paths for the includes and libraries



To restore the ibex library to its original state, simply overwrite the modified ibex_IntervalVector.h with the ibex_IntervalVector-original.h attached
