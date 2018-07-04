#Author Caleb Amoa Buahin
#Email caleb.buahin@gmail.com
#Date 2018
#License GNU General Public License (see <http: //www.gnu.org/licenses/> for details).
#The RHEComponent is a component-based radiative heat exchange stream transport model

TEMPLATE = lib
VERSION = 1.0.0
TARGET = RHEComponent
QT -= gui
QT += testlib


#Added for faster compilation
*msvc* { # visual studio spec filter
      QMAKE_CXXFLAGS += /MP /O2
  }


DEFINES += RHECOMPONENT_LIBRARY
DEFINES += USE_OPENMP
DEFINES += USE_MPI
DEFINES += USE_CVODE
DEFINES += USE_NETCDF
#DEFINES += USE_CVODE_OPENMP

#Compile as library or executable
contains(DEFINES,RHECOMPONENT_LIBRARY){
  TEMPLATE = lib
  message("Compiling RHEComponent as library")
} else {
  TEMPLATE = app
  CONFIG-=app_bundle
  message("Compiling RHEComponent as application")
}

CONFIG += c++11

linux{
CONFIG += debug_and_release
}

PRECOMPILED_HEADER = ./include/stdafx.h

INCLUDEPATH += .\
               ./include \
               ./../HydroCouple/include \
               ./../HydroCoupleSDK/include


HEADERS += ./include/stdafx.h\
           ./include/rhecomponent_global.h \
           ./include/rhecomponent.h \
           ./include/rhecomponentinfo.h \
           ./include/test/rhecomponenttest.h \
           ./include/elementinput.h \
           ./include/rhemodel.h \
           ./include/iboundarycondition.h \
           ./include/elementjunction.h \
           ./include/element.h \
           ./include/abstracttimeseriesbc.h \
           ./include/hydraulicstimeseriesbc.h \
           ./include/meteorologytimeseriesbc.h \
           ./include/radiativefluxtimeseriesbc.h \
           ./include/elementtimeseriesbc.h \
           ./include/elementoutput.h


SOURCES +=./src/stdafx.cpp \
          ./src/rhecomponent.cpp \
          ./src/rhecomponentinfo.cpp \
          ./src/main.cpp \
          ./src/elementinput.cpp \
          ./src/test/rhecomponenttest.cpp \
          ./src/rhemodel.cpp \
          ./src/elementjunction.cpp \
          ./src/element.cpp \
          ./src/element.cpp \
          ./src/rhecompute.cpp \
          ./src/rhemodelio.cpp \
          ./src/hydraulicstimeseriesbc.cpp \
          ./src/radiativefluxtimeseriesbc.cpp \
          ./src/abstracttimeseriesbc.cpp \
          ./src/elementtimeseriesbc.cpp \
          ./src/meteorologytimeseriesbc.cpp \
          ./src/elementoutput.cpp


macx{

    INCLUDEPATH += /usr/local \
                   /usr/local/include

    contains(DEFINES, USE_CVODE){
    message("CVODE enabled")
    LIBS += -L/usr/local/lib -lsundials_cvode
    }

    contains(DEFINES, USE_NETCDF){
    message("NetCDF enabled")
    LIBS += -L/usr/local/lib -lnetcdf-cxx4
    }

    contains(DEFINES,USE_OPENMP){

        QMAKE_CC = /usr/local/opt/llvm/bin/clang
        QMAKE_CXX = /usr/local/opt/llvm/bin/clang++
        QMAKE_LINK = /usr/local/opt/llvm/bin/clang++

        QMAKE_CFLAGS+= -fopenmp
        QMAKE_LFLAGS+= -fopenmp
        QMAKE_CXXFLAGS+= -fopenmp

        INCLUDEPATH += /usr/local/opt/llvm/lib/clang/5.0.0/include
        LIBS += -L /usr/local/opt/llvm/lib -lomp

      message("OpenMP enabled")
    } else {
      message("OpenMP disabled")
    }

    contains(DEFINES,USE_MPI){

        QMAKE_CXX = /usr/local/bin/mpicxx
        QMAKE_LINK = /usr/local/bin/mpicxx

        QMAKE_CFLAGS += $$system(mpicc --showme:compile)
        QMAKE_CXXFLAGS += $$system(mpic++ --showme:compile)
        QMAKE_LFLAGS += $$system(mpic++ --showme:link)

        LIBS += -L/usr/local/lib -lmpi

        message("MPI enabled")

    } else {

      message("MPI disabled")

    }
}

linux{

INCLUDEPATH += /usr/include \
               ../gdal/include

    contains(DEFINES,UTAH_CHPC){

         INCLUDEPATH += /uufs/chpc.utah.edu/sys/installdir/hdf5/1.8.17-c7/include \
                        /uufs/chpc.utah.edu/sys/installdir/netcdf-c/4.3.3.1/include \
                        /uufs/chpc.utah.edu/sys/installdir/netcdf-cxx/4.3.0-c7/include \
                        ../hypre/build/include


         LIBS += -L/uufs/chpc.utah.edu/sys/installdir/hdf5/1.8.17-c7/lib -lhdf5 \
                 -L/uufs/chpc.utah.edu/sys/installdir/netcdf-cxx/4.3.0-c7/lib -lnetcdf_c++4 \
                 -L../hypre/build/lib -lHYPRE

         message("Compiling on CHPC")
    }

    contains(DEFINES,USE_OPENMP){

    QMAKE_CFLAGS += -fopenmp
    QMAKE_LFLAGS += -fopenmp
    QMAKE_CXXFLAGS += -fopenmp

    LIBS += -L/usr/lib/x86_64-linux-gnu -lgomp

      message("OpenMP enabled")

    } else {

      message("OpenMP disabled")

    }
}

win32{

    #Windows vspkg package manager installation path
    VSPKGDIR = C:/vcpkg/installed/x64-windows

    INCLUDEPATH += $${VSPKGDIR}/include \
                   $${VSPKGDIR}/include/gdal

    CONFIG(debug, debug|release) {
    LIBS += -L$${VSPKGDIR}/debug/lib -lgdald
        } else {
    LIBS += -L$${VSPKGDIR}/lib -lgdal
    }

    contains(DEFINES, USE_CVODE){
    message("CVODE enabled")
    CONFIG(debug, debug|release) {
        LIBS += -L$${VSPKGDIR}/debug/lib -lsundials_cvode
        } else {
        LIBS += -L$${VSPKGDIR}/lib -lsundials_cvode
        }
    }

    contains(DEFINES, USE_NETCDF){
    message("NetCDF enabled")
    CONFIG(release, debug|release) {
        LIBS += -L$${VSPKGDIR}/lib -lnetcdf \
                -L$${VSPKGDIR}/lib -lnetcdf-cxx4
        } else {
        LIBS += -L$${VSPKGDIR}/debug/lib -lnetcdf \
                -L$${VSPKGDIR}/debug/lib -lnetcdf-cxx4
        }
    }

    contains(DEFINES,USE_OPENMP){

        QMAKE_CFLAGS += /openmp
        QMAKE_LFLAGS += /openmp
        QMAKE_CXXFLAGS += /openmp
        QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
        QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS

        message("OpenMP enabled")

     } else {

      message("OpenMP disabled")

     }

    contains(DEFINES,USE_MPI){

       LIBS += -L$$(MSMPI_LIB64)/ -lmsmpi

       message("MPI enabled")

     } else {

      message("MPI disabled")

     }
}

CONFIG(debug, debug|release) {

   DESTDIR = ./build/debug
   OBJECTS_DIR = $$DESTDIR/.obj
   MOC_DIR = $$DESTDIR/.moc
   RCC_DIR = $$DESTDIR/.qrc
   UI_DIR = $$DESTDIR/.ui

   macx{

    QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/build/debug/*HydroCoupleSDK.* ./build/debug/";
    LIBS += -L./../HydroCoupleSDK/build/debug -lHydroCoupleSDK.1.0.0

    }

   linux{

    QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/build/debug/*HydroCoupleSDK.* ./build/debug/";
    LIBS += -L./../HydroCoupleSDK/build/debug -lHydroCoupleSDK.so.1.0.0

    }

   win32{

    QMAKE_POST_LINK += "copy ./../HydroCoupleSDK/build/debug/*HydroCoupleSDK.* ./build/debug/";
    LIBS += -L./../HydroCoupleSDK/build/debug -lHydroCoupleSDK1

    }
}

CONFIG(release, debug|release) {

    RELEASE_EXTRAS = ./build/release
    OBJECTS_DIR = $$RELEASE_EXTRAS/.obj
    MOC_DIR = $$RELEASE_EXTRAS/.moc
    RCC_DIR = $$RELEASE_EXTRAS/.qrc
    UI_DIR = $$RELEASE_EXTRAS/.ui

   macx{
    LIBS += -L./../HydroCoupleSDK/lib/macx -lHydroCoupleSDK.1.0.0
    }

   linux{
    LIBS += -L./../HydroCoupleSDK/lib/linux -lHydroCoupleSDK.so.1.0.0
    }

   win32{
    LIBS += -L./../HydroCoupleSDK/lib/win32 -lHydroCoupleSDK1
    }

     contains(DEFINES,RHEComponent_LIBRARY){
         #MacOS
         macx{
             DESTDIR = lib/macx
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/macx/*HydroCoupleSDK.* ./lib/macx/";
        }

         #Linux
         linux{
             DESTDIR = lib/linux
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/linux/*HydroCoupleSDK.* ./lib/linux/";
        }

         #Windows
         win32{
             DESTDIR = lib/win32
             QMAKE_POST_LINK += "copy ./../HydroCoupleSDK/lib/win32/*HydroCoupleSDK* ./lib/win32/";
        }
    } else {
         #MacOS
         macx{
             DESTDIR = bin/macx
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/macx/*HydroCoupleSDK.* ./bin/macx/";
        }

         #Linux
         linux{
             DESTDIR = bin/linux
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/linux/*HydroCoupleSDK.* ./bin/linux/";
        }

         #Windows
         win32{
             DESTDIR = bin/win32
             QMAKE_POST_LINK += "copy ./../HydroCoupleSDK/lib/win32/*HydroCoupleSDK* ./bin/win32/";
        }
    }
}
