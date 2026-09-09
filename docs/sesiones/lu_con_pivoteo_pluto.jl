### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ e906c53a-8f53-44f7-a8b6-4f19f142186f
md"""
# Laboratorio: LU con pivoteo parcial

Primero construiremos **a mano** un ejemplo ``3\times 3`` para ver ``P``, ``L`` y ``U``. Después programaremos nuestro propio `mi_LU_pivoteo(A)` y resolveremos ``Ax=b``.

No usamos `LinearAlgebra` ni `A \ b`.
"""

# ╔═╡ 250c41c6-6311-4980-ad0d-21ea74d4db0b
md"""
## 1. Ejemplo ``3\times 3`` a mano

```math
A=\begin{pmatrix}
0&2&1\\
2&2&3\\
4&2&1
\end{pmatrix}.
```

El mayor elemento en valor absoluto de la primera columna es ``4``. Por pivoteo parcial intercambiamos las filas 1 y 3.
"""

# ╔═╡ 3bb5983d-056f-4bd2-ab5f-2be39b04da17
A = [
    0.0 2.0 1.0
    2.0 2.0 3.0
    4.0 2.0 1.0
]

# ╔═╡ c696811a-b917-4a6d-a1fe-59484ba70017
P = [
    0.0 0.0 1.0
    0.0 1.0 0.0
    1.0 0.0 0.0
]

# ╔═╡ 30cfc33d-6f74-475a-aead-f2f483edc9f3
P * A

# ╔═╡ 3d0de7c6-3acf-40e5-a0ac-65bcaa484427
md"""
Sobre ``PA`` hacemos eliminación. Los multiplicadores son

```math
m_{21}=\frac12,\qquad m_{32}=2.
```

Así obtenemos

```math
L=\begin{pmatrix}
1&0&0\\
1/2&1&0\\
0&2&1
\end{pmatrix},
\qquad
U=\begin{pmatrix}
4&2&1\\
0&1&5/2\\
0&0&-4
\end{pmatrix}.
```
"""

# ╔═╡ a0794cb6-cdbd-442c-829c-fb80adaa5afe
begin
    L_manual = [
        1.0 0.0 0.0
        0.5 1.0 0.0
        0.0 2.0 1.0
    ]
    U_manual = [
        4.0 2.0 1.0
        0.0 1.0 2.5
        0.0 0.0 -4.0
    ]
end

# ╔═╡ 74ca0a74-b671-4cad-8e53-84c03beadf68
md"""
## 2. La comprobación desnuda

Ejecuta ambas celdas y compara.
"""

# ╔═╡ 7a1270a0-fd09-4433-97af-413255a770dd
P * A

# ╔═╡ 67371a8d-e8e0-46bc-9668-6dba0d045b71
L_manual * U_manual

# ╔═╡ 8293e49d-9398-435d-baa5-f3284d92dabe
md"""
Debemos observar

```math
\boxed{PA=LU}.
```
"""

# ╔═╡ a486a8f5-10f7-4da0-be6e-06844065bd5d
md"""
## 3. Resolver ``Ax=b`` a mano

Tomemos

```math
x_{\rm exacta}=\begin{pmatrix}1\\2\\-1\end{pmatrix}
```

y construyamos ``b=Ax_{\rm exacta}``. Como ``PA=LU``, resolvemos

```math
Ly=Pb,\qquad Ux=y.
```
"""

# ╔═╡ aa17179c-2413-4431-a484-75ca4fd7e6a6
begin
    x_exacta = [1.0, 2.0, -1.0]
    b = A * x_exacta
    Pb = P * b
end

# ╔═╡ a46b2abd-9fd5-4af6-b37c-15b2ee04b812
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

# ╔═╡ 6c79c341-e6a2-4159-8504-6f3de4dc02ad
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

# ╔═╡ 3fbeb8fb-9d8b-44ca-baa6-d2c1912e6186
begin
    y_manual = forward_substitution(L_manual, Pb)
    x_manual = backward_substitution(U_manual, y_manual)
end

# ╔═╡ e15e33ed-f735-44bc-b10c-9da47f0d8c15
(x_exacta=x_exacta, x_obtenida=x_manual, error=x_manual-x_exacta)

# ╔═╡ 03a035dd-afd5-4705-b894-8ebf43105830
md"""
## 4. Ahora que lo haga nuestra máquina

En cada etapa, el algoritmo busca el pivote de mayor magnitud, intercambia filas, registra la permutación y continúa la eliminación.
"""

# ╔═╡ 7bb50a54-ba83-4230-a2d3-5d543350d640
function identidad_manual(n)
    I = zeros(Float64, n, n)
    for i in 1:n
        I[i,i] = 1.0
    end
    I
end

# ╔═╡ 1723a62b-756f-4b19-8fac-d7e170f38b23
function mi_LU_pivoteo(A)
    n, m = size(A)
    n == m || error("A debe ser cuadrada")

    U = Float64.(copy(A))
    L = identidad_manual(n)
    P = identidad_manual(n)

    for k in 1:n-1
        p = k
        mayor = abs(U[k,k])

        for i in k+1:n
            if abs(U[i,k]) > mayor
                mayor = abs(U[i,k])
                p = i
            end
        end

        mayor > 1e-14 || error("Matriz singular o pivote demasiado pequeño")

        if p != k
            U[[k,p], :] = U[[p,k], :]
            P[[k,p], :] = P[[p,k], :]
            if k > 1
                L[[k,p], 1:k-1] = L[[p,k], 1:k-1]
            end
        end

        for i in k+1:n
            m_ik = U[i,k] / U[k,k]
            L[i,k] = m_ik
            for j in k:n
                U[i,j] -= m_ik * U[k,j]
            end
        end
    end

    abs(U[n,n]) > 1e-14 || error("Matriz singular o pivote demasiado pequeño")
    P, L, U
end

# ╔═╡ a5cb69e0-0dca-4008-a87f-8e11ea19fad1
P_auto, L_auto, U_auto = mi_LU_pivoteo(A)

# ╔═╡ f1c44843-8858-4645-8a56-687d96dfcc84
(P=P_auto, L=L_auto, U=U_auto)

# ╔═╡ 36431190-8354-48d9-b4f7-b1d2a03146a0
begin
    izquierda = P_auto * A
    derecha = L_auto * U_auto
    (PA=izquierda, LU=derecha, diferencia=izquierda-derecha)
end

# ╔═╡ 5f8cc913-dd94-4913-b252-83d490ed83a6
md"""
## 5. Nuestro solucionador completo

```math
Ax=b\Longrightarrow Ly=Pb\Longrightarrow Ux=y.
```
"""

# ╔═╡ 8b95e653-4587-499a-893e-24f1ffbf641a
function resolver_LU_pivoteo(A, b)
    P, L, U = mi_LU_pivoteo(A)
    Pb = P * b
    y = forward_substitution(L, Pb)
    backward_substitution(U, y)
end

# ╔═╡ 32a09fa3-bf25-4f25-9eb0-8a6294148520
x_auto = resolver_LU_pivoteo(A, b)

# ╔═╡ 036e7b23-2da3-4309-999e-c61d7f3d568e
(x_exacta=x_exacta, x_obtenida=x_auto, error=x_auto-x_exacta)

# ╔═╡ e002efff-b9fe-4fda-9714-9d0ab0e4ac97
md"""
## 6. El pivote pequeño

Ahora probemos

```math
A_\varepsilon=
\begin{pmatrix}
\varepsilon&1\\
1&-1
\end{pmatrix},
\qquad \varepsilon=10^{-16}.
```

**Predice antes de ejecutar:** ¿qué fila escogerá nuestro pivoteo parcial?
"""

# ╔═╡ 9b2e7876-d7dd-4bac-ba7b-af2e9df72639
begin
    ε = 1e-16
    Aε = [ε 1.0; 1.0 -1.0]
    xε_exacta = [1.0, 1.0]
    bε = Aε * xε_exacta
end

# ╔═╡ 48b3311a-9255-4635-b40a-3a6057b253ca
Pε, Lε, Uε = mi_LU_pivoteo(Aε)

# ╔═╡ 32b64536-f99b-42d2-8a62-5d72254ae52d
(P=Pε, L=Lε, U=Uε)

# ╔═╡ d697eccc-5c67-423b-b007-e63a872a41ec
xε = resolver_LU_pivoteo(Aε, bε)

# ╔═╡ df64eb4c-ca87-4f5b-b7ea-329ae8884805
(x_exacta=xε_exacta, x_obtenida=xε, error=xε-xε_exacta)

# ╔═╡ 57b3d7c0-3d5a-4306-a4e8-39c9739bc1fe
md"""
## 7. Para discutir

1. ¿Qué información contiene ``P``?
2. ¿Por qué, en etapas posteriores, debemos cuidar lo que ya está almacenado en ``L``?
3. ¿Por qué el pivoteo parcial evita un multiplicador gigantesco cuando ``\varepsilon`` es muy pequeño?
4. ¿Cómo conecta este experimento con los errores de punto flotante?

> **Predice → ejecuta → observa → explica.**
"""

# ╔═╡ Cell order:
# ╟─e906c53a-8f53-44f7-a8b6-4f19f142186f
# ╟─250c41c6-6311-4980-ad0d-21ea74d4db0b
# ╠═3bb5983d-056f-4bd2-ab5f-2be39b04da17
# ╠═c696811a-b917-4a6d-a1fe-59484ba70017
# ╠═30cfc33d-6f74-475a-aead-f2f483edc9f3
# ╟─3d0de7c6-3acf-40e5-a0ac-65bcaa484427
# ╠═a0794cb6-cdbd-442c-829c-fb80adaa5afe
# ╟─74ca0a74-b671-4cad-8e53-84c03beadf68
# ╠═7a1270a0-fd09-4433-97af-413255a770dd
# ╠═67371a8d-e8e0-46bc-9668-6dba0d045b71
# ╟─8293e49d-9398-435d-baa5-f3284d92dabe
# ╟─a486a8f5-10f7-4da0-be6e-06844065bd5d
# ╠═aa17179c-2413-4431-a484-75ca4fd7e6a6
# ╠═a46b2abd-9fd5-4af6-b37c-15b2ee04b812
# ╠═6c79c341-e6a2-4159-8504-6f3de4dc02ad
# ╠═3fbeb8fb-9d8b-44ca-baa6-d2c1912e6186
# ╠═e15e33ed-f735-44bc-b10c-9da47f0d8c15
# ╟─03a035dd-afd5-4705-b894-8ebf43105830
# ╠═7bb50a54-ba83-4230-a2d3-5d543350d640
# ╠═1723a62b-756f-4b19-8fac-d7e170f38b23
# ╠═a5cb69e0-0dca-4008-a87f-8e11ea19fad1
# ╠═f1c44843-8858-4645-8a56-687d96dfcc84
# ╠═36431190-8354-48d9-b4f7-b1d2a03146a0
# ╟─5f8cc913-dd94-4913-b252-83d490ed83a6
# ╠═8b95e653-4587-499a-893e-24f1ffbf641a
# ╠═32a09fa3-bf25-4f25-9eb0-8a6294148520
# ╠═036e7b23-2da3-4309-999e-c61d7f3d568e
# ╟─e002efff-b9fe-4fda-9714-9d0ab0e4ac97
# ╠═9b2e7876-d7dd-4bac-ba7b-af2e9df72639
# ╠═48b3311a-9255-4635-b40a-3a6057b253ca
# ╠═32b64536-f99b-42d2-8a62-5d72254ae52d
# ╠═d697eccc-5c67-423b-b007-e63a872a41ec
# ╠═df64eb4c-ca87-4f5b-b7ea-329ae8884805
# ╟─57b3d7c0-3d5a-4306-a4e8-39c9739bc1fe
