project('gn', 'cpp','c',
        default_options : ['cpp_std=c++11',
			   'buildtype=debug',
			   'warning_level=3' ])

inc = include_directories('include')
mpi_dep_cpp = dependency('mpi', language: 'cpp')
mpi_dep_c = dependency('mpi', language: 'c')

umfpack = meson.get_compiler('cpp').find_library('umfpack',
  required: true)

# umfpack_dep = dependency('umfpack')
executable('basic', ['main.cc',
          					'src/utilities.cc',
                    'src/sparseDiracMatrixClass.cc',
                    'src/importanceSamplingIntegrator.cc',
                    'src/gaussianRandomGenerator.cc',
                    'src/determinantWrapper.cc',
                    'src/ranlxd.c',
                    'src/suiteSparseDet.c'
                    ], 
           include_directories: inc,
           dependencies: [mpi_dep_c,mpi_dep_cpp,umfpack])
