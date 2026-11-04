module SVAnimExt
    using StructuralVibration, Statistics, Makie, Meshes

    import StructuralVibration: animate_mesh, viz_mesh, fast_cmap

    include("mesh/common_mesh.jl")
    include("anim/anim_utils.jl")
end