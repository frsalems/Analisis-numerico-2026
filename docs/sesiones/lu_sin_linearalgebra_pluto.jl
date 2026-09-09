### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ c51c1ca5-1afb-45a6-a3c1-e17721daf415
md"""
# Resolver \(Ax=b\) con nuestro propio LU

Aquí no usamos `LinearAlgebra` ni `A\b`.

Construiremos nosotros mismos

```math
A=LU,
```

y resolveremos

```math
Ly=b,\qquad Ux=y.
```

> **Meta:** ver el algoritmo por dentro: factorizar → sustituir → verificar.
"""

# ╔═╡ 5a7967b4-b40a-43f6-8426-328c97b8ff34
md"""
## 1. Sustitución hacia adelante

Para una matriz triangular inferior ``L``, resolvemos ``Ly=b`` de arriba hacia abajo.
"""

# ╔═╡ 8f8bcf97-5e5b-46dd-ac2c-a742f91e0290
function forward_substitution(L, b)
    n = length(b)
    y = zeros(Float64, n)
    for i in 1:n
        s = 0.0
        for j in 1:i-1
            s += L[i,j] * y[j]
        end
        y[i] = (b[i] - s) / L[i,i]
    end
    y
end

# ╔═╡ fb77151f-47cb-4dcc-9c57-c76fda24df47
md"""
## 2. Sustitución hacia atrás

Para una matriz triangular superior ``U``, resolvemos ``Ux=y`` de abajo hacia arriba.
"""

# ╔═╡ 0a9e539b-1858-4bce-8bbd-ed92e48196ba
function backward_substitution(U, y)
    n = length(y)
    x = zeros(Float64, n)
    for i in n:-1:1
        s = 0.0
        for j in i+1:n
            s += U[i,j] * x[j]
        end
        x[i] = (y[i] - s) / U[i,i]
    end
    x
end

# ╔═╡ e1bff2a5-8413-44ec-a0b2-919ab18bf95c
function identidad_manual(n)
    L = zeros(Float64, n, n)
    for i in 1:n
        L[i,i] = 1.0
    end
    L
end

# ╔═╡ 5f9166a6-1a76-4ffb-9707-91711e7a455e
md"""
## 3. Nuestra factorización LU

Usamos eliminación gaussiana **sin pivoteo**. Los multiplicadores

```math
m_{ik}=\frac{U_{ik}}{U_{kk}}
```

se guardan debajo de la diagonal de ``L``.
"""

# ╔═╡ ba97e391-919b-4e08-8036-fc4a59fe791e
function mi_LU(A)
    n, m = size(A)
    n == m || error("A debe ser cuadrada.")

    U = Float64.(copy(A))
    L = identidad_manual(n)

    for k in 1:n-1
        abs(U[k,k]) > 1e-14 ||
            error("Pivote cero o muy pequeño: esta versión no usa pivoteo.")

        for i in k+1:n
            mult = U[i,k] / U[k,k]
            L[i,k] = mult
            for j in k:n
                U[i,j] -= mult * U[k,j]
            end
        end
    end

    abs(U[n,n]) > 1e-14 ||
        error("Pivote cero o muy pequeño: esta versión no usa pivoteo.")

    L, U
end

# ╔═╡ f75887b0-c7c3-457b-8dbf-dfab0c68ef48
md"""
## 4. El solucionador completo

La función sigue exactamente la cadena

```math
A\to(L,U)\to y\to x.
```
"""

# ╔═╡ 96816118-eace-46b1-a499-903fc7d7fec3
function resolver_LU(A, b)
    n, m = size(A)
    n == m || error("A debe ser cuadrada.")
    length(b) == n || error("Dimensiones incompatibles.")

    L, U = mi_LU(A)
    y = forward_substitution(L, Float64.(b))
    x = backward_substitution(U, y)

    x, L, U, y
end

# ╔═╡ 53d39858-e6a8-405e-b367-4c6a2c9f7008
md"""
## 5. Ejemplo (`4 X 4)`

Usaremos la misma matriz de nuestro ejemplo de factorización:

```math
A=
\begin{pmatrix}
2&1&-1&2\\
4&5&0&3\\
-2&8&11&-4\\
2&-5&3&11
\end{pmatrix}.
```

Escogemos un vector ``b`` y dejamos que nuestro algoritmo haga todo el trabajo.
"""

# ╔═╡ efcd2966-ea89-4e48-abe6-f7c15ac6afc7
begin 
A = [
     2.0   1.0  -1.0   2.0
     4.0   5.0   0.0   3.0
    -2.0   8.0  11.0  -4.0
     2.0  -5.0   3.0  11.0
]

b = [4.0, 12.0, 7.0, 20.0]
end	

# ╔═╡ 451a877d-2a4e-4cf9-b1bc-64af300069b0
x, L, U, y = resolver_LU(A, b)

# ╔═╡ a181f3e8-d359-475e-9108-dbc485632332
md"""### La matriz ``L``"""

# ╔═╡ ac655773-cff5-496c-b64d-6ffe1a54bd78
L

# ╔═╡ 0f1c4193-42a5-4df5-b837-3baad9f2f844
md"""### La matriz ``U``"""

# ╔═╡ f7a1f22b-06c8-4b29-8bdc-87de8d19a81c
U

# ╔═╡ 03b269fa-8b3e-4ff3-9a34-7d76ea497054
md"""### Primero resolvemos ``Ly=b`` y obtenemos"""

# ╔═╡ 98812f56-d440-41e0-af04-e199437207af
y

# ╔═╡ c7a24aab-8ff2-45e6-b3ed-32f6f04b1866
md"""### Después resolvemos ``Ux=y`` y obtenemos"""

# ╔═╡ b17b4430-91f0-4cab-a53d-ebd77e3a3fd8
x

# ╔═╡ b60b6df3-379d-47ad-bcd1-c90aaec59a75
md"""
## 6. Verificación sin `LinearAlgebra`

Calculamos nosotros mismos

```math
r=b-Ax.
```
"""

# ╔═╡ 4ed2170d-2cef-424e-b20a-837f2708c9f4
function matvec_manual(A, x)
    n, m = size(A)
    length(x) == m || error("Dimensiones incompatibles.")
    z = zeros(Float64, n)
    for i in 1:n
        for j in 1:m
            z[i] += A[i,j] * x[j]
        end
    end
    z
end

# ╔═╡ db819dca-6d51-44ef-a499-c9b9f302980b
residuo = b - matvec_manual(A, x)

# ╔═╡ 853aea20-f04b-4b36-9846-1aaf367c5d8b
md"""
Si los componentes del residuo son cero o muy pequeños, la solución calculada
satisface numéricamente el sistema.

> **Ojo:** residuo pequeño y error pequeño no son la misma afirmación.
"""

# ╔═╡ ba1555e3-c44e-408d-9a3d-df1c2a9bedfb
md"""
## 7. Cuando nuestro LU se topa con pared

Consideremos ahora el sistema

```math
\begin{pmatrix}
0 & 1\\
1 & 1
\end{pmatrix}
\begin{pmatrix}
x_1\\
x_2
\end{pmatrix}
=
\begin{pmatrix}
1\\
2
\end{pmatrix}.

```
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.10"
manifest_format = "2.0"
project_hash = "da39a3ee5e6b4b0d3255bfef95601890afd80709"

[deps]
"""

# ╔═╡ Cell order:
# ╠═c51c1ca5-1afb-45a6-a3c1-e17721daf415
# ╠═5a7967b4-b40a-43f6-8426-328c97b8ff34
# ╠═8f8bcf97-5e5b-46dd-ac2c-a742f91e0290
# ╠═fb77151f-47cb-4dcc-9c57-c76fda24df47
# ╠═0a9e539b-1858-4bce-8bbd-ed92e48196ba
# ╠═e1bff2a5-8413-44ec-a0b2-919ab18bf95c
# ╠═5f9166a6-1a76-4ffb-9707-91711e7a455e
# ╠═ba97e391-919b-4e08-8036-fc4a59fe791e
# ╠═f75887b0-c7c3-457b-8dbf-dfab0c68ef48
# ╠═96816118-eace-46b1-a499-903fc7d7fec3
# ╠═53d39858-e6a8-405e-b367-4c6a2c9f7008
# ╠═efcd2966-ea89-4e48-abe6-f7c15ac6afc7
# ╠═451a877d-2a4e-4cf9-b1bc-64af300069b0
# ╟─a181f3e8-d359-475e-9108-dbc485632332
# ╠═ac655773-cff5-496c-b64d-6ffe1a54bd78
# ╠═0f1c4193-42a5-4df5-b837-3baad9f2f844
# ╠═f7a1f22b-06c8-4b29-8bdc-87de8d19a81c
# ╠═03b269fa-8b3e-4ff3-9a34-7d76ea497054
# ╠═98812f56-d440-41e0-af04-e199437207af
# ╠═c7a24aab-8ff2-45e6-b3ed-32f6f04b1866
# ╠═b17b4430-91f0-4cab-a53d-ebd77e3a3fd8
# ╠═b60b6df3-379d-47ad-bcd1-c90aaec59a75
# ╠═4ed2170d-2cef-424e-b20a-837f2708c9f4
# ╠═db819dca-6d51-44ef-a499-c9b9f302980b
# ╠═853aea20-f04b-4b36-9846-1aaf367c5d8b
# ╠═ba1555e3-c44e-408d-9a3d-df1c2a9bedfb
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
