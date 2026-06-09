import numpy as np
import matplotlib.pyplot as plt

# Parâmetros do mecanismo (cm)
L1 = 10.0
L2 = 13.0
L3 = 8.0
L4 = 10.0

# Parâmetros do método
TOL_REL = 1e-6
MAX_ITER = 100
NUM_BETA = 600
BETA_MIN = 0.0
BETA_MAX = 2 * np.pi

# Constante da equação (independe de x e β)
K = (L1**2 + L2**2 - L3**2 + L4**2) / (2 * L2 * L4)


def f(x, beta):
    return (L1/L2)*np.cos(beta) - (L1/L4)*np.cos(x) - np.cos(beta - x) + K


def df(x, beta):
    return (L1/L4)*np.sin(x) - np.sin(beta - x)


def newton_raphson(x0, beta, tol=TOL_REL, max_iter=MAX_ITER):
    x_n = x0
    for i in range(1, max_iter + 1):
        fx = f(x_n, beta)
        dfx = df(x_n, beta)

        # if abs(dfx) < 1e-14:
        #     return x_n, i, False

        x_new = x_n - fx / (max(abs(dfx), 1e-14) * np.sign(dfx))

        erro = abs((x_new - x_n) / x_new) if abs(x_new) > 1e-14 else abs(x_new - x_n)
        x_n = x_new

        if erro < tol:
            return x_n, i, True

    return x_n, max_iter, False


def resolver(betas, x0):
    n = len(betas)
    solucoes = np.full(n, np.nan)
    iteracoes = np.zeros(n, dtype=int)
    convergencias = np.zeros(n, dtype=bool)

    for i, beta in enumerate(betas):
        x_sol, n_iter, conv = newton_raphson(x0, beta)

        iteracoes[i] = n_iter
        convergencias[i] = conv

        if conv:
            solucoes[i] = np.mod(x_sol, 2*np.pi)
            # x0 = solucoes[i]

    return solucoes, iteracoes, convergencias


def raizes_por_varredura(beta, sementes):
    # Acha as raízes distintas em um dado β testando várias sementes fixas.
    # Usado só para "ancorar" o início de cada ramo (não usa forma fechada).
    raizes = []
    for s in sementes:
        x, _, conv = newton_raphson(s, beta)
        if conv:
            xm = np.mod(x, 2*np.pi)
            nova = all(min(abs(xm - r), 2*np.pi - abs(xm - r)) > 1e-3 for r in raizes)
            if nova:
                raizes.append(xm)
    return sorted(raizes)


def tracar_ramos(betas):
    # Rastreia os dois ramos α(β) com Newton-Raphson de SEMENTE MÓVEL (continuação):
    # cada β é resolvido usando como semente a solução convergida do β anterior.
    # Isso mantém o método "preso" a um ramo e evita o salto caótico entre raízes.
    # Quando entra na zona proibida o Newton falha (NaN); ao reaparecer solução,
    # re-ancora com uma varredura de sementes.
    n = len(betas)
    ramo1 = np.full(n, np.nan)
    ramo2 = np.full(n, np.nan)
    sementes = np.linspace(0, 2*np.pi, 16, endpoint=False)

    g1 = None   # semente corrente do ramo 1 (valor anterior)
    g2 = None   # semente corrente do ramo 2

    for i, beta in enumerate(betas):
        if g1 is None or g2 is None:
            raizes = raizes_por_varredura(beta, sementes)
            
            
            if len(raizes) >= 2:
                g1, g2 = raizes[0], raizes[-1]
            elif len(raizes) == 1:
                g1 = g2 = raizes[0]
            else:
                continue   # ainda na zona proibida

        x1, _, c1 = newton_raphson(g1, beta)
        x2, _, c2 = newton_raphson(g2, beta)

        if c1:
            ramo1[i] = np.mod(x1, 2*np.pi)
            g1 = ramo1[i]
        else:
            g1 = None

        if c2:
            ramo2[i] = np.mod(x2, 2*np.pi)
            g2 = ramo2[i]
        else:
            g2 = None

    return ramo1, ramo2


def plotar(betas,
           solucoes1, iteracoes1, convergencias1,
           solucoes2, iteracoes2, convergencias2):

    betas_deg = np.degrees(betas)

    solucoes1_deg = np.degrees(solucoes1)
    solucoes2_deg = np.degrees(solucoes2)

    fig, ax = plt.subplots(figsize=(10, 6))

    fig.suptitle(
        "Equação de Freudenstein — Newton-Raphson",
        fontsize=14,
        fontweight="bold"
    )

    # Regiões onde pelo menos uma das inicializações falhou
    falhas = ~(convergencias1 & convergencias2)

    inicio = None

    for i in range(len(falhas)):
        if falhas[i] and inicio is None:
            inicio = i

        terminou = (
            inicio is not None and
            (
                not falhas[i] or
                i == len(falhas) - 1
            )
        )

        if terminou:
            fim = i if not falhas[i] else i + 1

            ax.axvspan(
                betas_deg[inicio],
                betas_deg[fim - 1],
                color="red",
                alpha=0.15
            )

            inicio = None

    ax.plot(
        betas_deg,
        solucoes1_deg,
        "o-",
        markersize=4,
        linewidth=1.5,
        label=r"$x_0=-0.1$"
    )

    ax.plot(
        betas_deg,
        solucoes2_deg,
        "s-",
        markersize=4,
        linewidth=1.5,
        label=r"$x_0=2\pi/3$"
    )

    ax.set_xlabel(r"$\beta$ (graus)", fontsize=12)
    ax.set_ylabel(r"$x$ (graus)", fontsize=12)
    ax.set_title("Raízes obtidas para dois valores iniciais", fontsize=12)

    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.6)

    plt.tight_layout()
    plt.savefig(
        "freudenstein_resultados.png",
        dpi=150,
        bbox_inches="tight"
    )
    plt.show()

    

    print("\n✔ Gráfico salvo em 'freudenstein_resultados.png'")


def plotar_bacias(x0_1=-0.1, x0_2=2*np.pi/3):
    # Varre o plano (β, x0): para cada semente x0 e cada β aplica Newton-Raphson
    # e colore o ponto conforme o ramo encontrado. A fronteira entre as bacias
    # é fractal — é isso que torna a semente fixa x0=2π/3 instável em parte do
    # intervalo (a reta horizontal atravessa a região marmorizada).

    NB = 6400   # resolução em β
    NX = 4800   # resolução em x0

    betas = np.linspace(BETA_MIN, BETA_MAX, NB)
    x0s = np.linspace(-1.5, 6.5, NX)

    BB, XX = np.meshgrid(betas, x0s)

    # Newton-Raphson vetorizado (mesma iteração de newton_raphson, sem o laço escalar)
    x = XX.copy()
    with np.errstate(all="ignore"):
        for _ in range(MAX_ITER):
            x = x - f(x, BB) / df(x, BB)
        residuo = np.abs(f(x, BB))

    convergiu = residuo < 1e-6   # tolerância no resíduo |f|
    x_mod = np.mod(x, 2*np.pi)

    # Valor real de cada ramo via Newton de semente móvel (continuação),
    # usado para rotular em qual ramo o método caiu. Sem forma fechada.
    raiz1, raiz2 = tracar_ramos(betas)

    # Distância circular de cada solução às duas raízes (para escolher o ramo)
    d1 = np.abs(x_mod - raiz1[None, :]) % (2*np.pi)
    d1 = np.minimum(d1, 2*np.pi - d1)
    d2 = np.abs(x_mod - raiz2[None, :]) % (2*np.pi)
    d2 = np.minimum(d2, 2*np.pi - d2)
    perto_de_1 = d1 <= d2

    # Código de cor: 0 = não converge, 1 = ramo 1, 2 = ramo 2
    codigo = np.zeros(x_mod.shape, dtype=int)
    codigo[convergiu & perto_de_1] = 1
    codigo[convergiu & ~perto_de_1] = 2

    cores = {
        0: (0.07, 0.07, 0.07),   # preto   — não converge / sem raiz
        1: (0.12, 0.47, 0.71),   # azul    — ramo 1
        2: (1.00, 0.50, 0.05),   # laranja — ramo 2
    }
    imagem = np.zeros(codigo.shape + (3,))
    for valor, cor in cores.items():
        imagem[codigo == valor] = cor

    fig, ax = plt.subplots(figsize=(11, 6.5))

    fig.suptitle(
        "Equação de Freudenstein — Bacias de atração do Newton-Raphson",
        fontsize=14,
        fontweight="bold"
    )

    ax.imshow(
        imagem,
        origin="lower",
        aspect="auto",
        extent=[np.degrees(BETA_MIN), np.degrees(BETA_MAX), x0s.min(), x0s.max()],
        interpolation="nearest"
    )

    ax.axhline(
        x0_1,
        color="cyan",
        linewidth=1.8,
        linestyle="--",
        label=r"$x_0=-0.1$"
    )

    ax.axhline(
        x0_2,
        color="magenta",
        linewidth=1.8,
        linestyle="--",
        label=r"$x_0=2\pi/3$"
    )

    ax.set_xlabel(r"$\beta$ (graus)", fontsize=12)
    ax.set_ylabel(r"$x_0$ (semente, rad)", fontsize=12)
    ax.set_title("Preto = não converge · azul/laranja = ramo encontrado", fontsize=12)

    ax.legend(loc="upper right", framealpha=0.9)

    plt.tight_layout()
    plt.savefig(
        "freudenstein_bacias.png",
        dpi=150,
        bbox_inches="tight"
    )
    plt.show()

    from PIL import Image
    img_uint8 = (imagem * 255).astype(np.uint8)
    Image.fromarray(img_uint8).save("freudenstein_bacias_hires.png")

    print("\n✔ Gráfico salvo em 'freudenstein_bacias.png'")


def main():
    print(f"l1={L1}, l2={L2}, l3={L3}, l4={L4} cm")
    print(f"β ∈ [0°, {np.degrees(BETA_MAX):.1f}°] | {NUM_BETA} valores | tol={TOL_REL:.0e} | max_iter={MAX_ITER}")

    betas = np.linspace(BETA_MIN, BETA_MAX, NUM_BETA)
    solucoes1, iteracoes1, convergencias1 = resolver(betas, -0.1)
    solucoes2, iteracoes2, convergencias2 = resolver(betas, 2*np.pi/3)

    # # Tabela de resultados
    # print(f"\n{'β (°)':>10} {'x (°)':>10} {'x (rad)':>12} {'Iter':>6} {'Conv':>5}")
    # print("-" * 50)
    # for b, s, it, c in zip(betas, solucoes, iteracoes, convergencias):
    #     print(f"{np.degrees(b):10.3f} {np.degrees(s):10.3f} {s:12.6f} {it:6d} {'✔' if c else '✘':>5}")

    print(f"\nConvergiram (x0=-0.1): {np.sum(convergencias1)}/{len(betas)}")
    print(f"Convergiram (x0=2π/3): {np.sum(convergencias2)}/{len(betas)}")

    plotar(
        betas,
        solucoes1, iteracoes1, convergencias1,
        solucoes2, iteracoes2, convergencias2
    )

    plotar_bacias()


if __name__ == "__main__":
    main()