
const MBXjson   = joinpath(@__DIR__, "mbx.json")
const libmbx    = dlopen(joinpath(@__DIR__, "libmbx.so"), Libdl.RTLD_GLOBAL) 
const E_MBX2JMD = 0.043 # kcal/mol --> eV

# /**
#  * Initializes the system in the heap.
#  * @param[in] coords Pointer to the coordinates (size 3N)
#  * @param[in] nat_monomers Pointer to an array of the number of atoms in each monomer
#  * @param[in] at_names Pointer to an array with the atom names of all the whole system
#  * @param[in] monomers Pointer to the list of monomer ids in your system
#  * @param[in] nmon Number of monomers
#  * @param[in] json_file Name of the json configuration file
#  */
# void initialize_system_py_(double* coords, int* nat_monomers, char** at_names, char** monomers, int* nmon,
#                            char* json_file)

function mbx_initialize_system(vars)
  sym = dlsym(libmbx, :initialize_system_py_)

  @ccall $sym(
    vars.xyz::Ptr{Cdouble},
    vars.num_ats::Ptr{Cint},
    vars.at_nams::Ptr{Ptr{Cchar}},
    vars.mon_nams::Ptr{Ptr{Cchar}},
    vars.num_mon::Ref{Cint},
    vars.json::Ptr{Cchar}
  )::Cvoid
end

function mbx_initialize_system(xyz, num_ats, at_nams, mon_nams, num_mon, json)
  sym = dlsym(libmbx, :initialize_system_py_)

  @ccall $sym(
    xyz::Ptr{Cdouble},
    num_ats::Ptr{Cint},
    at_nams::Ptr{Ptr{Cchar}},
    mon_nams::Ptr{Ptr{Cchar}},
    num_mon::Ref{Cint},
    json::Ptr{Cchar}
  )::Cvoid
end

function mbx_get_energy(xyz, nats)
  sym = dlsym(libmbx, :get_energy_)

  E   = Ref{Cdouble}(0)

  @ccall $sym(xyz::Ptr{Cdouble}, nats::Ref{Cint}, E::Ref{Cdouble})::Cvoid

  E[] * E_MBX2JMD
end

function mbx_get_energy_grad!(G, xyz, nats)
  sym = dlsym(libmbx, :get_energy_g_)

  E   = Ref{Cdouble}(0)

  @ccall $sym(
    xyz::Ptr{Cdouble}, 
    nats::Ref{Cint}, 
    E::Ref{Cdouble},
    G::Ptr{Cdouble}
  )::Cvoid

  G .*= E_MBX2JMD

  E[] * E_MBX2JMD
end

function mbx_finalize_system()
  sym = dlsym(libmbx, :finalize_system_)

  @ccall $sym()::Cvoid
end