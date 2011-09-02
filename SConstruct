#
# builds static and shared library 'giraffe'
# and applications 'xcca' and 'testsuite'
#





##### define include directories
#includedirs = Split("""/reg/neh/home/feldkamp/hdf5/include
#  /reg/neh/home/feldkamp/fftw/include
#  /reg/neh/home/feldkamp/tiff/include
#  edf_support""")

include_dirs = Split("""
  /usr/local/hdf5/include
  /usr/local/include
  edf_support
  """)


lib_dirs = Split("""
  /usr/local/hdf5/lib
  /usr/local/lib
  edf_support
  """)

libs = Split("""
  hdf5
  fftw3
  tiff
  """)



##### define sources to be compiled
sources_giraffe = Split("""
  crosscorrelator.cpp
  arrayclasses.cpp
  arraydataIO.cpp
  fouriertransformer.cpp
  """)
  
sources_edf_support = Split("""
  edf_support/edf.cpp
  edf_support/edfkeywords.cpp
  edf_support/timer.cpp
  edf_support/util.cpp
  """)

sources_xcca = Split("""
  xcca.cpp
  analyze.cpp
  """)

sources_testsuite = Split("""
  testsuite.cpp
  tester.cpp
  """)






##### go to work
env = Environment()
env.Append( CPPPATH=include_dirs, LIBPATH=lib_dirs, LIBS=libs )

env.SharedLibrary('libgiraffe', sources_giraffe+sources_edf_support)
env.StaticLibrary('libgiraffe', sources_giraffe+sources_edf_support)

env.Program('xcca', sources_xcca+sources_giraffe+sources_edf_support)
env.Program('testsuite', sources_testsuite+sources_giraffe+sources_edf_support)
