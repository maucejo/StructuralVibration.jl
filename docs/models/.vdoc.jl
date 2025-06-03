#
#
#
#
#
#
#
#
#
#
#
#
#
#| echo: false
#| output: false
using StructuralVibration
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| echo: false
@doc Bar
```
#
#
#
#
#
#| echo: false
@doc Rod
```
#
#
#
#
#
#| echo: false
@doc Strings
```
#
#
#
#
#
#| echo: false
@doc Beam
```
#
#
#
#
#
#
#
#| echo: false
@doc modefreq(b::Bar)
```
#
#
#
#
#
#| echo: false
@doc modeshape(b::Bar)
```
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| echo: false
@doc Membrane
```
#
#
#
#
#
#| echo: false
@doc Plate
```
#
#
#
#
#
#
#
#| echo: false
@doc modefreq(p::Plate)
```
#
#
#
#
#
#| echo: false
@doc modeshape(p::Plate)
```
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| echo: false
@doc Sdof
```
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| echo: false
@doc Mdof
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| echo: false
@doc MdofMesh
```
#
#
#
#
#
#
#
#| echo: false

@doc assembly(model::Mdof)
```
#
#
#
#
#
#| echo: false

@doc apply_bc
```
#
#
#
#
#
#
#| echo: false

@doc eigenmode(K::Matrix{Float64}, M::Matrix{Float64})
```
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| echo: false

@doc OneDMesh
```
#
#
#
#
#
#
#
#
#| echo: false

@doc assembly(model::Beam, mesh::OneDMesh)
```
#
#
#
#
#
#
#| echo: false

@doc rayleigh_damping_matrix(K, M, α::Float64, β::Float64)
```
#
#
#
#
#
#
#| echo: false

@doc modal_damping_matrix
```
#
#
#
#
#
#
#| echo: false

@doc selection_matrix
```
#
#
#
#
#
#
#
#
#
#
#
