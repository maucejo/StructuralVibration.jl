
"""
    @show_struct struct StructName
        field1::Type1
        field2::Type2
        ...
    end

Macro that automatically generates an enhanced display method for a structure.

This macro applies to a structure definition and overloads the `Base.show` function
to provide formatted output that includes the structure name and, for each field,
its name, type, and value.

# Display Format
- For simple types (Number, String, Symbol): displays the value directly
- For arrays (AbstractArray):
  - If elements are simple types and the array contains 10 elements or less:
    displays all elements
  - If elements are simple types and the array contains more than 10 elements:
    displays the first 5 elements followed by "..."
  - Always includes an array summary (dimensions and type)
- For other types: displays the type and standard representation

# Example
```julia
@show_struct struct Person
    name::String
    age::Int
    scores::Vector{Float64}
end

p = Person("John", 30, [15.5, 17.0, 16.5])
println(p)
# Output:
#   Person
#     name: String John
#     age: Int64 30
#     scores: Vector{Float64} [15.5, 17.0, 16.5]
```
"""
macro show_struct(expr)
    if expr.head != :struct
        error("@show_struct macro must be used with a struct definition")
    end

    # Extract struct information
    struct_def = expr.args[2]

    # Extract the base type name without parameters
    base_name = nothing
    if isa(struct_def, Symbol)
        # Simple case: struct Name
        base_name = struct_def
    elseif isa(struct_def, Expr)
        if struct_def.head == :<:
            # Case: struct Name <: Parent or struct Name{T} <: Parent
            subtype_expr = struct_def.args[1]
            if isa(subtype_expr, Symbol)
                base_name = subtype_expr
            elseif isa(subtype_expr, Expr) && subtype_expr.head == :curly
                # Handle parametric type: struct Name{T}
                base_name = subtype_expr.args[1]
            end
        elseif struct_def.head == :curly
            # Case: struct Name{T}
            base_name = struct_def.args[1]
        end
    end

    if base_name === nothing
        error("Unsupported struct definition format")
    end

    # Définition normale de la structure
    result = :($(esc(expr)))

    # Une seule méthode générique pour gérer tous les cas
    show_def = quote
        function Base.show(io::IO, obj::T) where {T <: $(esc(base_name))}
            println(io, "  $(typeof(obj))")

            for field in fieldnames(typeof(obj))
                value = getfield(obj, field)
                value_type = typeof(value)

                if value_type <: Number || value_type <: AbstractString || value_type <: Symbol
                    println(io, "  "^(2) * "$field: $value_type $value")
                elseif value_type <: AbstractArray
                    val_str = "..."
                    if !isempty(value) && (eltype(value) <: Number || eltype(value) <: AbstractString || eltype(value) <: Symbol)
                        if length(value) <= 10
                            val_str = "$value"
                        else
                            val_str = string(value[1:5]) * "..."
                        end
                    end
                    println(io, "  "^(2) * "$field: $(summary(value)) $val_str")
                elseif isstructtype(value_type)
                    # Pour les structures, afficher uniquement le type, pas la valeur
                    println(io, "  "^(2) * "$field: $value_type")
                else
                    println(io, "  "^(2) * "$field: $value_type $value")
                end
            end
        end
    end

    # Combine the structure definition and the show method
    return quote
        $result
        $show_def
    end
end