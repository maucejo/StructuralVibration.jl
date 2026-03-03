module SVMeshExt
    using StructuralVibration, Meshes, Statistics

    import StructuralVibration: detrend_mesh, renumber_element_connectivity,
        build_mesh, deformed_grid

    include("common_mesh.jl")
    include("mesh_utils.jl")
end