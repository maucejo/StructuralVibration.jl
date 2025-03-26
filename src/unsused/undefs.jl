# Convenient function to create arrays of undefs
"""
    undefs(::Type{T}, dims::Int...) where T
    undefs(dims::Int...)

Create an array of undefs of type `T` and dimensions `dims`. If `T` is not provided, the default type is `Float64`.

# Inputs
- `T::Type`: Type of the array
- `dims::Int...`: Dimensions of the array

# Output
- `arr`: Array of undefs
"""
undefs(::Type{T}, dims::Int...) where T = Array{T}(undef, dims)
undefs(::Type{T}, dims::Tuple{Vararg{Int}}) where T = Array{T}(undef, dims...)
undefs(dims::Int...) = Array{Float64}(undef, dims)
undefs(dims::Tuple{Vararg{Int}}) = Array{Float64}(undef, dims)