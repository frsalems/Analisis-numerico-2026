### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ c7e83cc0-4e31-4786-a342-b09c407f7e86
begin
    using LinearAlgebra
end

# ╔═╡ 2c6d1fa2-8f7d-4cf6-8c54-27dbf9f2c301
md"""
# Laboratorio — QR a mano

La idea de este cuaderno es **no llamar a `qr(A)` hasta el final**.

Vamos a construir

```math
A=QR
```

mediante Gram-Schmidt, verificarla y usarla para resolver un problema de mínimos cuadrados.

La computadora será calculadora y verificadora. La estructura matemática la ponemos nosotros.
"""

# ╔═╡ 4fdc5bb2-c9de-43fc-b3a2-1923b375ce79
md"""
## 1. La matriz

Trabajaremos con

```math
A=\begin{bmatrix}
1&1\\
1&0\\
0&1
\end{bmatrix}.
```

Antes de ejecutar nada:

1. ¿son ortogonales sus columnas?
2. ¿tienen norma uno?
3. ¿qué tendremos que hacer para convertirlas en una base ortonormal?
"""

# ╔═╡ 2c1dfa0c-f62f-4f34-a0ac-7d68ea39a1e8
A = [1.0 1.0;
     1.0 0.0;
     0.0 1.0]

# ╔═╡ 13189d8c-ec48-4394-a2e4-ac2fc1d5ce02
begin
    a1 = A[:,1]
    a2 = A[:,2]
    (dot(a1,a2), norm(a1), norm(a2))
end

# ╔═╡ f62cbd17-6302-4ea7-a0a0-e2952aaf8c3b
md"""
## 2. Primera dirección

```math
q_1=\frac{a_1}{\|a_1\|},
\qquad
r_{11}=\|a_1\|.
```
"""

# ╔═╡ 33529924-8776-44b1-b14c-78854687a1ca
begin
    r11 = norm(a1)
    q1 = a1/r11
    q1
end

# ╔═╡ fbf4a04b-fbf1-491e-a2a8-39bef0be6edb
dot(q1,q1)

# ╔═╡ 52e31aef-b2cb-4a4e-9f2a-1b7392482d40
md"""
## 3. Quitamos a $a_2$ su componente en dirección $q_1$

```math
r_{12}=q_1^Ta_2,
```

```math
v_2=a_2-r_{12}q_1.
```

Antes de ejecutar: ¿por qué $v_2$ debe ser perpendicular a $q_1$?
"""

# ╔═╡ edd3ae4b-4c73-4bfa-8fb0-56accfe8dd19
begin
    r12 = dot(q1,a2)
    v2 = a2-r12*q1
    (r12,v2,dot(q1,v2))
end

# ╔═╡ 56d75006-2d0d-443c-ac06-ad408b59dd80
md"""
## 4. Segunda dirección

```math
r_{22}=\|v_2\|,
\qquad
q_2=\frac{v_2}{r_{22}}.
```
"""

# ╔═╡ b2e7afbb-7346-414d-852d-57d381a79550
begin
    r22 = norm(v2)
    q2 = v2/r22
    (r22,q2)
end

# ╔═╡ f237371f-1677-4f1b-8211-d328579c0c54
(dot(q1,q2),dot(q2,q2))

# ╔═╡ ac98f894-5e19-4bdc-9891-105601107c97
md"""
## 5. Construimos $Q$ y $R$

```math
Q=\begin{bmatrix}q_1&q_2\end{bmatrix},
\qquad
R=\begin{bmatrix}r_{11}&r_{12}\\0&r_{22}\end{bmatrix}.
```
"""

# ╔═╡ e8d1c1c3-fe3c-47b9-95e8-fac78cf19ee0
begin
    Q = hcat(q1,q2)
    R = [r11 r12;
         0.0 r22]
    (Q,R)
end

# ╔═╡ 2362bb33-9752-403f-95fb-067f85b2f6e0
md"""
## 6. Verificación

Debemos tener

```math
Q^TQ\approx I,
\qquad
QR\approx A.
```
"""

# ╔═╡ 8152bdc3-e95b-4776-b89f-50a2ca759576
begin
    error_ortonormalidad = norm(Q'*Q-I)
    error_reconstruccion = norm(A-Q*R)
    (error_ortonormalidad,error_reconstruccion)
end

# ╔═╡ 2686c5ce-d25b-462f-a2ef-2c60892d90d8
md"""
## 7. Un problema rectangular

Tomemos

```math
b=\begin{bmatrix}2\\1\\1\end{bmatrix}.
```

Buscamos el $x$ que minimiza

```math
\|b-Ax\|_2.
```

Como $A=QR$, resolveremos

```math
Rx=Q^Tb.
```
"""

# ╔═╡ 6b331b9a-8b9b-4817-b0ed-b22c20a7d121
b = [2.0,1.0,1.0]

# ╔═╡ 1beef06d-e57d-452d-9383-5b8c06405ac8
z = Q'*b

# ╔═╡ e54d8d48-d612-4e38-aaad-46dc7e88c0b0
md"""
Como $R$ es triangular superior,

```math
x_2=\frac{z_2}{r_{22}},
```

```math
x_1=\frac{z_1-r_{12}x_2}{r_{11}}.
```
"""

# ╔═╡ 04cfbf10-588f-4669-944a-88d4b67e12a5
begin
    x2 = z[2]/R[2,2]
    x1 = (z[1]-R[1,2]*x2)/R[1,1]
    x_manual = [x1,x2]
end

# ╔═╡ 9de6d751-f1b3-4607-8c17-c9ad2f36b31d
md"""
## 8. El residuo

Debe cumplirse

```math
A^Tr\approx0,
\qquad
r=b-Ax.
```
"""

# ╔═╡ e061c7c6-3c89-4598-917e-51c2b8bb7d11
begin
    r = b-A*x_manual
    (x_manual,r,A'*r,norm(r))
end

# ╔═╡ bf2c4d66-06c8-4c38-96a7-d67f029246bd
md"""
## 9. Ahora sí: dejamos entrar a Julia

Sólo ahora comparamos con `qr(A)` y `A \\ b`.
"""

# ╔═╡ 040fa0de-1d8e-47de-a41a-10cd767b8cbb
begin
    F = qr(A)

    Qj_full = Matrix(F.Q)
    Qj = Qj_full[:,1:size(A,2)]
    Rj = Matrix(F.R)

    x_julia = A \ b

    println("Nuestra solución = ", x_manual)
    println("A\\b              = ", x_julia)

    println("||A-Qj*Rj|| = ", norm(A-Qj*Rj))
    println("||QQ'-QjQj'|| = ", norm(Q*Q' - Qj*Qj'))
    println("diferencia soluciones = ", norm(x_manual-x_julia))
end

# ╔═╡ ebd83ad8-9068-480d-a3ca-b585648c46bf
md"""
## 10. Reto 1

Trabaja con

```math
B=\begin{bmatrix}
1&1\\
1&2\\
1&3\\
1&4
\end{bmatrix}.
```

Sin llamar a `qr(B)`:

1. construye $q_1$;
2. elimina de la segunda columna su proyección sobre $q_1$;
3. construye $q_2$;
4. forma $Q$ y $R$;
5. verifica $Q^TQ\approx I$;
6. verifica $QR\approx B$.
"""

# ╔═╡ f43607bb-7122-43f4-8df6-acde1527e106
B = [1.0 1.0;
     1.0 2.0;
     1.0 3.0;
     1.0 4.0]

# ╔═╡ 16bef331-6376-4923-9d26-13164311661f
md"""
## 11. Reto 2 — QR y ajuste de una recta

Toma

```julia
t = [0.,1.,2.,3.,4.]
y = [1.0,1.8,3.2,3.9,5.1]
```

Construye la matriz de diseño $X$ y después:

1. factoriza $X=QR$ manualmente;
2. calcula $Q^Ty$;
3. resuelve $R\beta=Q^Ty$;
4. calcula los residuos;
5. comprueba $X^Tr\approx0$;
6. compara al final con `X \\ y`.

**No uses `qr(X)` hasta haber terminado.**
"""

# ╔═╡ 15b47dc7-943d-4d5a-9f41-25df3b24064b
md"""
## Cierre

Cada paso de Gram-Schmidt hace una cosa geométrica concreta:

```math
\text{quitar proyecciones}
\longrightarrow
\text{crear direcciones ortogonales}
\longrightarrow
\text{construir }Q
\longrightarrow
\text{obtener }R.
```

Y después:

```math
\text{mínimos cuadrados}
\longrightarrow
Rx=Q^Tb.
```

La próxima pregunta será numérica:

> ¿Es Gram-Schmidt siempre una buena manera de calcular QR?

Ahí empieza la historia de Householder.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.10"
manifest_format = "2.0"
project_hash = "ac1187e548c6ab173ac57d4e72da1620216bce54"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╟─2c6d1fa2-8f7d-4cf6-8c54-27dbf9f2c301
# ╠═c7e83cc0-4e31-4786-a342-b09c407f7e86
# ╠═4fdc5bb2-c9de-43fc-b3a2-1923b375ce79
# ╠═2c1dfa0c-f62f-4f34-a0ac-7d68ea39a1e8
# ╠═13189d8c-ec48-4394-a2e4-ac2fc1d5ce02
# ╟─f62cbd17-6302-4ea7-a0a0-e2952aaf8c3b
# ╠═33529924-8776-44b1-b14c-78854687a1ca
# ╠═fbf4a04b-fbf1-491e-a2a8-39bef0be6edb
# ╟─52e31aef-b2cb-4a4e-9f2a-1b7392482d40
# ╠═edd3ae4b-4c73-4bfa-8fb0-56accfe8dd19
# ╟─56d75006-2d0d-443c-ac06-ad408b59dd80
# ╠═b2e7afbb-7346-414d-852d-57d381a79550
# ╠═f237371f-1677-4f1b-8211-d328579c0c54
# ╟─ac98f894-5e19-4bdc-9891-105601107c97
# ╠═e8d1c1c3-fe3c-47b9-95e8-fac78cf19ee0
# ╟─2362bb33-9752-403f-95fb-067f85b2f6e0
# ╠═8152bdc3-e95b-4776-b89f-50a2ca759576
# ╟─2686c5ce-d25b-462f-a2ef-2c60892d90d8
# ╠═6b331b9a-8b9b-4817-b0ed-b22c20a7d121
# ╠═1beef06d-e57d-452d-9383-5b8c06405ac8
# ╟─e54d8d48-d612-4e38-aaad-46dc7e88c0b0
# ╠═04cfbf10-588f-4669-944a-88d4b67e12a5
# ╟─9de6d751-f1b3-4607-8c17-c9ad2f36b31d
# ╠═e061c7c6-3c89-4598-917e-51c2b8bb7d11
# ╟─bf2c4d66-06c8-4c38-96a7-d67f029246bd
# ╠═040fa0de-1d8e-47de-a41a-10cd767b8cbb
# ╟─ebd83ad8-9068-480d-a3ca-b585648c46bf
# ╠═f43607bb-7122-43f4-8df6-acde1527e106
# ╟─16bef331-6376-4923-9d26-13164311661f
# ╟─15b47dc7-943d-4d5a-9f41-25df3b24064b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
