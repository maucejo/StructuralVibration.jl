module SVAnimExt
    using StructuralVibration, Statistics, Makie, Meshes

    import StructuralVibration: animate_mesh, viz_mesh, fast_cmap

    include("common_mesh.jl")
    include("anim_utils.jl")
end