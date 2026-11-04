module SVMeshExt
    using StructuralVibration, Meshes, Statistics

    import StructuralVibration: detrend_mesh, renumber_element_connectivity,
        build_mesh, deformed_grid

    include("mesh/common_mesh.jl")
    include("mesh/mesh_utils.jl")
end