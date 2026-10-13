### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 72778e10-e01f-431f-9689-26dd31986042
begin
    using Markdown
    using InteractiveUtils
    using PlutoUI
    using Plots
end

# ╔═╡ 566a432c-9e75-45a0-9c0c-e6799b3b5851
md"""
# Métodos no lineales: una mirada geométrica

Este cuaderno acompaña la **Sesión 08**.

La idea no es comenzar con algoritmos, sino **ver qué están haciendo**.

Primero estudiaremos la geometría de la bisección y después veremos
cómo el método de Newton utiliza rectas tangentes para construir
aproximaciones sucesivas a una raíz.
"""

# ╔═╡ e70b37b9-e877-4a1b-926f-d4c694884872
TableOfContents()

# ╔═╡ 5972d5ae-a8c7-42c2-8600-2c3f0a91d955
begin
    f(x) = x^3 - x - 2
    df(x) = 3x^2 - 1
end

# ╔═╡ efb3dc7a-d0a1-40c3-9e14-d9f097eb732c
md"""
# 1. Bisección

Queremos resolver

```math
f(x)=x^3-x-2=0.
```

Tomamos inicialmente

```math
a_0=1,\qquad b_0=2.
```

Como

```math
f(1)<0,\qquad f(2)>0,
```

la función cambia de signo en el intervalo $[1,2]$.

Por continuidad existe al menos una raíz entre ambos puntos.

La bisección consiste en observar repetidamente el punto medio

```math
c_k=\frac{a_k+b_k}{2}
```

y conservar la mitad del intervalo en la que continúa existiendo
un cambio de signo.
"""

# ╔═╡ 9989c541-1f7f-48ee-9680-18fccbce740c
function historia_biseccion(f, a, b, n)
    historia = NamedTuple[]

    for k in 0:n
        c = (a + b)/2
        fa, fb, fc = f(a), f(b), f(c)

        push!(historia, (
            k = k,
            a = a,
            b = b,
            c = c,
            fa = fa,
            fb = fb,
            fc = fc
        ))

        k == n && break

        if fa*fc < 0
            b = c
        else
            a = c
        end
    end

    historia
end

# ╔═╡ 75a9a15e-3ec2-4a02-bb30-dd07a9944f2b
@bind k_bis Slider(0:10, default=0, show_value=true)

# ╔═╡ a0884fa7-e4b4-4bf4-9088-262a83384cda
begin
    hist_bis = historia_biseccion(f, 1.0, 2.0, k_bis)
    paso_bis = hist_bis[end]

    a = paso_bis.a
    b = paso_bis.b
    c = paso_bis.c
    fa = paso_bis.fa
    fb = paso_bis.fb
    fc = paso_bis.fc
end

# ╔═╡ 7eaeb702-b8d6-4564-82a0-1a61f630cd9a
md"""
## Iteración $k_bis$

En esta iteración,

```math
a=$(round(a,digits=8)),
\qquad
b=$(round(b,digits=8)),
\qquad
c=\frac{a+b}{2}=$(round(c,digits=8)).
```

Además,

```math
f(a)=$(round(fa,digits=8)),
\qquad
f(b)=$(round(fb,digits=8)),
\qquad
f(c)=$(round(fc,digits=8)).
```

El producto de los valores en los extremos es

```math
f(a)f(b)=$(round(fa*fb,digits=8)).
```

Como es negativo, los extremos están a lados distintos del eje horizontal.

$(
    fa*fc < 0 ?
    "El cambio de signo está entre a y c, así que el siguiente intervalo será [a,c]." :
    "El cambio de signo está entre c y b, así que el siguiente intervalo será [c,b]."
)
"""

# ╔═╡ 38286261-9ff0-43fd-8d41-5d17276b071b
begin
    xs = range(0.5, 2.2, length=600)

    p_bis = plot(
        xs,
        f.(xs),
        linewidth=2,
        xlabel="x",
        ylabel="f(x)",
        label="f(x)=x³-x-2",
        title="Bisección: el cambio de signo encierra la raíz"
    )

    hline!(p_bis, [0], linestyle=:dash, label="y=0")

    scatter!(p_bis, [a], [fa], markersize=7, label="a, f(a)")
    scatter!(p_bis, [b], [fb], markersize=7, label="b, f(b)")
    scatter!(p_bis, [c], [fc], markersize=7, label="c, f(c)")

    plot!(p_bis, [a,a], [0,fa], linestyle=:dot, label="")
    plot!(p_bis, [b,b], [0,fb], linestyle=:dot, label="")
    plot!(p_bis, [c,c], [0,fc], linestyle=:dot, label="")

    p_bis
end

# ╔═╡ 2f5bfb5c-f4c5-4535-b5a8-dcc927edc3fc
begin
    hist10 = historia_biseccion(f, 1.0, 2.0, 10)
    ks_bis = [h.k for h in hist10]
    anchos = [h.b - h.a for h in hist10]

    plot(
        ks_bis,
        anchos,
        marker=:circle,
        linewidth=2,
        xlabel="iteración",
        ylabel="longitud del intervalo",
        label="b-a",
        title="Bisección: el intervalo se reduce a la mitad"
    )
end

# ╔═╡ 64818a27-8e6a-45fe-9740-948e3dd16284
md"""
# 2. Newton: la tangente encuentra la siguiente aproximación

Supongamos que $x_k$ es una aproximación a la raíz.

La recta tangente a $f$ en $x_k$ es

```math
y=f(x_k)+f'(x_k)(x-x_k).
```

La siguiente aproximación $x_{k+1}$ se obtiene haciendo $y=0$:

```math
0=f(x_k)+f'(x_k)(x_{k+1}-x_k).
```

Por tanto,

```math
\boxed{
x_{k+1}
=
x_k-\frac{f(x_k)}{f'(x_k)}
}
```

La construcción es geométrica: **la tangente en $x_k$ corta al eje horizontal
en $x_{k+1}$**.
"""

# ╔═╡ b0357429-e8b8-461a-a2c3-abd802221cad
@bind x0_newton Slider(1.1:0.05:2.0, default=1.5, show_value=true)

# ╔═╡ 3b41e02b-3982-4eb9-a2fe-6c28af5fbcdd
@bind n_newton Slider(0:6, default=1, show_value=true)

# ╔═╡ 916edc22-25a0-49fb-b386-4e70f898df14
function historia_newton(f, df, x0, n; derivada_min=1e-12)
    xs = [float(x0)]

    for _ in 1:n
        x = xs[end]
        dfx = df(x)

        abs(dfx) < derivada_min && break

        push!(xs, x - f(x)/dfx)
    end

    xs
end

# ╔═╡ 3596410d-34a9-49d5-b141-807716773643
xs_newton = historia_newton(f, df, x0_newton, n_newton)

# ╔═╡ a31e85c7-6b63-43f2-a69c-f23519ef7fb7
begin
    dominio = range(0.5, 2.5, length=800)

    p_newton = plot(
        dominio,
        f.(dominio),
        linewidth=2,
        xlabel="x",
        ylabel="f(x)",
        label="f(x)=x³-x-2",
        title="Newton: las tangentes construyen las aproximaciones"
    )

    hline!(p_newton, [0], linestyle=:dash, label="y=0")

    for j in 1:(length(xs_newton)-1)
        xk = xs_newton[j]
        xnext = xs_newton[j+1]
        fxk = f(xk)

        plot!(
            p_newton,
            [xk, xnext],
            [fxk, 0],
            linewidth=2,
            label=(j == 1 ? "tangentes" : "")
        )

        plot!(
            p_newton,
            [xk, xk],
            [0, fxk],
            linestyle=:dot,
            label=""
        )

        scatter!(p_newton, [xk], [fxk], markersize=6, label="")
        scatter!(p_newton, [xnext], [0], markersize=6, label="")
    end

    p_newton
end

# ╔═╡ c94e6a50-db9f-47eb-9285-060bfbca6350
begin
    filas_newton = [
        (
            iteración = k-1,
            x = x,
            f_x = f(x),
            residuo = abs(f(x))
        )
        for (k,x) in enumerate(xs_newton)
    ]

    filas_newton
end

# ╔═╡ cabf2bb3-30c2-40cc-9ad6-b612efeda23f
md"""
## ¿Qué conviene observar?

Mueva $x_0$ y aumente poco a poco el número de iteraciones.

Observe que cada paso contiene dos movimientos geométricos:

1. desde el eje horizontal subimos hasta la gráfica;
2. desde ese punto seguimos la tangente hasta volver a cortar el eje horizontal.

Cuando $x_0$ está razonablemente cerca de la raíz, las aproximaciones
se acercan muy rápidamente.

Pero esto todavía no cuenta toda la historia.

Newton también puede:

- entrar en ciclos;
- dar saltos enormes;
- fallar cuando $f'(x_k)$ es muy pequeña.

Esos comportamientos serán el tema de un segundo cuaderno:
**Newton y sus patologías**.
"""

# ╔═╡ 64804e8a-aedc-404f-afeb-de36058197bb
begin
    raiz_ref = 1.5213797068045676

    hist_b = historia_biseccion(f, 1.0, 2.0, 12)
    aprox_b = [h.c for h in hist_b]
    err_b = abs.(aprox_b .- raiz_ref)

    xs_n = historia_newton(f, df, 1.5, 6)
    err_n = abs.(xs_n .- raiz_ref)

    p_comp = plot(
        0:(length(err_b)-1),
        max.(err_b, eps()),
        linewidth=2,
        marker=:circle,
        yscale=:log10,
        xlabel="iteración",
        ylabel="error absoluto",
        label="bisección",
        title="Bisección y Newton: rapidez de convergencia"
    )

    plot!(
        p_comp,
        0:(length(err_n)-1),
        max.(err_n, eps()),
        linewidth=2,
        marker=:circle,
        label="Newton"
    )

    p_comp
end

# ╔═╡ 27102b39-dd5f-4b0c-944a-2e2e43a96fcd
md"""
# Para pensar

1. ¿Qué información utiliza bisección que Newton no necesita?
2. ¿Por qué bisección es lenta pero robusta?
3. ¿Qué representa geométricamente la fórmula de Newton?
4. ¿Por qué una derivada muy pequeña puede producir un salto enorme?
5. ¿Qué ventaja tendría combinar la seguridad de bisección con la rapidez de Newton?
"""

# ╔═╡ Cell order:
# ╠═72778e10-e01f-431f-9689-26dd31986042
# ╟─566a432c-9e75-45a0-9c0c-e6799b3b5851
# ╠═e70b37b9-e877-4a1b-926f-d4c694884872
# ╠═5972d5ae-a8c7-42c2-8600-2c3f0a91d955
# ╟─efb3dc7a-d0a1-40c3-9e14-d9f097eb732c
# ╠═9989c541-1f7f-48ee-9680-18fccbce740c
# ╠═75a9a15e-3ec2-4a02-bb30-dd07a9944f2b
# ╠═a0884fa7-e4b4-4bf4-9088-262a83384cda
# ╟─7eaeb702-b8d6-4564-82a0-1a61f630cd9a
# ╠═38286261-9ff0-43fd-8d41-5d17276b071b
# ╠═2f5bfb5c-f4c5-4535-b5a8-dcc927edc3fc
# ╟─64818a27-8e6a-45fe-9740-948e3dd16284
# ╠═b0357429-e8b8-461a-a2c3-abd802221cad
# ╠═3b41e02b-3982-4eb9-a2fe-6c28af5fbcdd
# ╠═916edc22-25a0-49fb-b386-4e70f898df14
# ╠═3596410d-34a9-49d5-b141-807716773643
# ╠═a31e85c7-6b63-43f2-a69c-f23519ef7fb7
# ╠═c94e6a50-db9f-47eb-9285-060bfbca6350
# ╟─cabf2bb3-30c2-40cc-9ad6-b612efeda23f
# ╠═64804e8a-aedc-404f-afeb-de36058197bb
# ╟─27102b39-dd5f-4b0c-944a-2e2e43a96fcd
