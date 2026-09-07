### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 88fed3e5-ebe9-4beb-ae1b-52a9f5a33a89
begin
    using LinearAlgebra
    using PlutoUI
end

# ╔═╡ 4c520319-2000-4351-a617-db1518e3d572
md"""
# ⚡ Laboratorio interactivo: circuito de tres mallas

Aquí seguiremos la cadena

```math
\text{circuito} \longrightarrow \text{leyes físicas} \longrightarrow AI=b
\longrightarrow \text{solución numérica}.
```

La regla del laboratorio es: **predice → mueve → observa → explica**.
"""

# ╔═╡ 6de125e0-e101-4702-9a4c-5095df9b1cc5
md"""
## 1. Recordatorio físico

**Ley de Ohm**

```math
V=RI.
```

**Primera ley de Kirchhoff:** en un nodo,

```math
\sum I=0.
```

**Segunda ley de Kirchhoff:** en una malla cerrada,

```math
\sum V=0.
```

Usaremos la ley de voltajes de Kirchhoff junto con la ley de Ohm.
"""

# ╔═╡ 08fcef25-d5c1-4afb-beb1-c12d88d08e05
md"""
## 2. Circuito de tres mallas

Sean ``I_1``, ``I_2`` e ``I_3`` las corrientes de malla, orientadas en sentido horario.

Tomaremos

```math
R_1=4\,\Omega,\quad R_2=6\,\Omega,\quad R_3=5\,\Omega,
```

```math
R_{12}=2\,\Omega,\qquad R_{23}=3\,\Omega.
```

Las fuentes iniciales serán ``E_1=12V``, ``E_2=5V`` y ``E_3=10V``.
"""

# ╔═╡ ec19a1f5-b189-4151-a54b-a8c449cd978c
LocalResource("circuito_tres_mallas.png")

# ╔═╡ 9c8612bb-85f4-4d8a-aed5-368b678d4226
md"""
## 3. De Kirchhoff al sistema lineal

```math
(R_1+R_{12})I_1-R_{12}I_2=E_1,
```

```math
-R_{12}I_1+(R_2+R_{12}+R_{23})I_2-R_{23}I_3=E_2,
```

```math
-R_{23}I_2+(R_3+R_{23})I_3=E_3.
```

Por tanto,

```math
\begin{pmatrix}
6&-2&0\\
-2&11&-3\\
0&-3&8
\end{pmatrix}
\begin{pmatrix}I_1\\I_2\\I_3\end{pmatrix}
=
\begin{pmatrix}E_1\\E_2\\E_3\end{pmatrix}.
```
"""

# ╔═╡ cddc7411-ee8c-41d7-821a-7841bdc17ab5
md"""
## 4. Mueve las fuentes ⚡

**``E_1``:** $(@bind E1 Slider(0.0:0.5:20.0, default=12.0, show_value=true)) V

**``E_2``:** $(@bind E2 Slider(0.0:0.5:20.0, default=5.0, show_value=true)) V

**``E_3``:** $(@bind E3 Slider(0.0:0.5:20.0, default=10.0, show_value=true)) V
"""

# ╔═╡ e92c1125-1037-42b3-8d38-22af77d0d1f8
begin
    R1, R2, R3 = 4.0, 6.0, 5.0
    R12, R23 = 2.0, 3.0
    A = [
        R1 + R12   -R12              0.0
        -R12       R2 + R12 + R23   -R23
         0.0       -R23              R3 + R23
    ]
    b = [E1, E2, E3]
    I = A \ b
end

# ╔═╡ 3c03c9ad-614e-42d2-a6d0-2ee144f9f347

md"""
### Corrientes

- **I₁ = $(round(I[1], digits=4)) A**
- **I₂ = $(round(I[2], digits=4)) A**
- **I₃ = $(round(I[3], digits=4)) A**

Mueve una sola fuente y observa cómo responden las tres corrientes.
"""

# ╔═╡ 2b9572e8-7b57-441d-ad81-88cc9d21e94a
begin
    residuo_relativo = norm(A*I-b) / max(norm(b), eps())
    numero_condicion = cond(A)
    (; residuo_relativo, numero_condicion)
end

# ╔═╡ 6137e180-021a-433d-b2d4-9743cb593621
md"""
### ¿Qué nos dicen estos números?

- **Residuo relativo:** $(round(residuo_relativo, sigdigits=4))
- **Número de condición:** $(round(numero_condicion, digits=4))

El residuo relativo es extremadamente pequeño: la solución calculada satisface
el sistema con gran precisión numérica.

El número de condición es moderado, por lo que este sistema no presenta
la sensibilidad extrema que observamos en las matrices de Hilbert.

> **Observa:** resolver ``A I = b`` no termina al obtener ``I``.
> También debemos preguntarnos qué tan confiable y sensible es la solución.
"""

# ╔═╡ 1b6dd957-995c-4392-8245-e1b9701996a2
md"""
## 5. La estructura de ``A``

``A`` es simétrica. Los elementos diagonales reúnen resistencias de cada malla y los términos fuera de la diagonal representan resistencias compartidas.

La matriz contiene, por tanto, información sobre la **topología del circuito**.
"""

# ╔═╡ bea9a20f-f4e2-4209-912a-6b53742fcc60
(; simetrica = issymmetric(A), autovalores = eigvals(Symmetric(A)))

# ╔═╡ 9cfa5e30-1175-4f16-b753-7682a2cde224
md"""
## 6. Experimento

1. Mantén ``E_2`` y ``E_3`` fijas y aumenta ``E_1``.
2. **Predice** primero cómo cambiarán ``I_1``, ``I_2`` e ``I_3``.
3. Mueve la regleta y compara.
4. Busca valores para los cuales alguna corriente cambie de signo.
5. Explica físicamente ese cambio.

> **Pregunta final:** ¿qué parte del modelo está contenida en ``A`` y cuál en ``b``?
"""

# ╔═╡ 77b9d71d-6589-4acb-8b1b-7a39484b4583
md"""
## 7. Siguiente observacion 

Aquí las resistencias permanecen fijas. Más adelante podemos convertir también las resistencias en regletas.

Entonces cambiará la propia matriz ``A`` y podremos observar simultáneamente cómo cambian las corrientes y el número de condición.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.10"
manifest_format = "2.0"
project_hash = "640bbf5877ef1f5e9e4848ce7b0954f33a8a7931"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

    [deps.ColorTypes.weakdeps]
    StyledStrings = "f489334b-da3d-4c2e-b8f0-e476e12c162b"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"
"""

# ╔═╡ Cell order:
# ╟─4c520319-2000-4351-a617-db1518e3d572
# ╠═88fed3e5-ebe9-4beb-ae1b-52a9f5a33a89
# ╟─6de125e0-e101-4702-9a4c-5095df9b1cc5
# ╟─08fcef25-d5c1-4afb-beb1-c12d88d08e05
# ╠═ec19a1f5-b189-4151-a54b-a8c449cd978c
# ╟─9c8612bb-85f4-4d8a-aed5-368b678d4226
# ╟─cddc7411-ee8c-41d7-821a-7841bdc17ab5
# ╠═e92c1125-1037-42b3-8d38-22af77d0d1f8
# ╠═3c03c9ad-614e-42d2-a6d0-2ee144f9f347
# ╠═2b9572e8-7b57-441d-ad81-88cc9d21e94a
# ╠═6137e180-021a-433d-b2d4-9743cb593621
# ╟─1b6dd957-995c-4392-8245-e1b9701996a2
# ╠═bea9a20f-f4e2-4209-912a-6b53742fcc60
# ╟─9cfa5e30-1175-4f16-b753-7682a2cde224
# ╟─77b9d71d-6589-4acb-8b1b-7a39484b4583
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
