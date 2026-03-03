module SVAnimExt
    using StructuralVibration, LinearAlgebra, Makie, Meshes

    import StructuralVibration: animate_mesh, viz_mesh

    include("common_mesh.jl")
    include("anim_utils.jl")
end