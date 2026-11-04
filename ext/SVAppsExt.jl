module SVAppsExt
    using StructuralVibration, LinearAlgebra, GLMakie

    import StructuralVibration: free_response_1dof, forced_response_1dof, free_response_2dof, forced_response_2dof, newmark_2dof

    include("apps/apps_utils.jl")
    include("apps/apps_1dof.jl")
    include("apps/apps_2dof.jl")
end