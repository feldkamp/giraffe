Help("""
-----------------------------------------------------------------------------
 builds static and shared library 'giraffe'
 and applications 'xcca' and 'testsuite' that apply giraffe functionality
 
-----------------------------------------------------------------------------
 options: (the first in the list is the default)
    host=psexport    :   SLAC's psexport
    host=jan         :   Jan's computer
    
    debug=0          :   no debug flags, optimize code
    debug=1          :   compile with debug flags --> slower code
    
    static=true      :   build static library 'giraffe.a'
    shared=true      :   build shared library 'giraffe.so'
    xcca=true        :   build executable 'xcca'
    testsuite=true   :   build executable 'testsuite'
-----------------------------------------------------------------------------
""")

host=ARGUMENTS.get('host', 'psexport')
debug=ARGUMENTS.get('debug', 0)


# define what to build
build_giraffe_static_lib=ARGUMENTS.get('static', True)
build_giraffe_shared_lib=ARGUMENTS.get('shared', True)
build_xcca=ARGUMENTS.get('xcca', True)
build_testsuite=ARGUMENTS.get('testsuite', True)
build_filehandler=ARGUMENTS.get('filehandler', True)




##### define include directories
if (host=='psexport'):
   include_dirs = Split("""
      /reg/neh/home/feldkamp/fftw/include
      /reg/neh/home/feldkamp/hdf5/include
      /reg/neh/home/feldkamp/tiff/include
      /reg/neh/home/feldkamp/packages_c/boost_1_47_0
      edf_support
      """)
   lib_dirs = Split("""
      /reg/neh/home/feldkamp/fftw/lib
      /reg/neh/home/feldkamp/hdf5/lib
      /reg/neh/home/feldkamp/tiff/lib
      /reg/neh/home/feldkamp/packages_c/boost_1_47_0/stage/lib
      .
      """)
elif (host=='jan'):
   include_dirs = Split("""
      /usr/local/hdf5/include
      /usr/local/include
      edf_support
      """)
   lib_dirs = Split("""
      /usr/local/hdf5/lib
      /usr/local/lib
      .
      """)
else:
   print "No valid host has been defined (host=", host, ")"
   print "Some important include and libraries paths may be missing"




##### define external libraries to be included
libs = Split("""
   hdf5
   fftw3
   tiff
   libboost_program_options
   """)


##### define sources to be compiled
sources_giraffe = Split("""
   crosscorrelator.cpp
   arrayclasses.cpp
   arraydataIO.cpp
   fouriertransformer.cpp
   """)
  
sources_edf = Split("""
   edf_support/edf.cpp
   edf_support/edfkeywords.cpp
   edf_support/timer.cpp
   edf_support/util.cpp
   """)

sources_xcca = Split("""
  xcca.cpp
  """)

sources_testsuite = Split("""
  testsuite.cpp
  tester.cpp
  """)

sources_filehandler = Split("""
  filehandler.cpp
  """)



##### go to work
print "using host >>>", host, "<<< (change this with option host=<value>)"

env = Environment(CPPPATH=include_dirs, LIBPATH=lib_dirs, LIBS=libs )


if debug:
   env.Append(CCFLAGS = "-g")
else:
   env.Append(CCFLAGS = "-O3")


objs_giraffe_static = env.StaticObject(sources_giraffe)
objs_giraffe_shared = env.SharedObject(sources_giraffe)
objs_edf_static = env.StaticObject(sources_edf)
objs_edf_shared = env.SharedObject(sources_edf)
objs_xcca = env.StaticObject(sources_xcca)
objs_testsuite = env.StaticObject(sources_testsuite)
objs_filehandler = env.StaticObject(sources_filehandler)


env.StaticLibrary( 'libgiraffe_static', objs_giraffe_static+objs_edf_static )
env.SharedLibrary( 'libgiraffe_shared', objs_giraffe_shared+objs_edf_shared, RPATH=lib_dirs )

env.Program( 'xcca', objs_xcca+objs_giraffe_static+objs_edf_static )
env.Program( 'testsuite', objs_testsuite+objs_giraffe_static+objs_edf_static )
env.Program( 'filehandler', objs_filehandler+objs_giraffe_static+objs_edf_static,RPATH=lib_dirs )
