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
NUM_BETA = 60
BETA_MIN = 0.0
BETA_MAX = 5 * np.pi / 6

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

        if abs(dfx) < 1e-14:
            return x_n, i, False

        x_new = x_n - fx / dfx

        erro = abs((x_new - x_n) / x_new) if abs(x_new) > 1e-14 else abs(x_new - x_n)
        x_n = x_new

        if erro < tol:
            return x_n, i, True

    return x_n, max_iter, False


def resolver(betas):
    n = len(betas)
    solucoes = np.zeros(n)
    iteracoes = np.zeros(n, dtype=int)
    convergencias = np.zeros(n, dtype=bool)

    x0 = betas[0] if betas[0] != 0 else 0.5

    for i, beta in enumerate(betas):
        x_sol, n_iter, conv = newton_raphson(x0, beta)
        solucoes[i], iteracoes[i], convergencias[i] = x_sol, n_iter, conv
        x0 = x_sol  # continuação numérica

    return solucoes, iteracoes, convergencias


def plotar(betas, solucoes, iteracoes):
    betas_deg = np.degrees(betas)
    solucoes_deg = np.degrees(solucoes)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("Equação de Freudenstein — Newton-Raphson", fontsize=14, fontweight="bold")

    ax1.plot(betas_deg, solucoes_deg, "o-", color="#2563EB", markersize=4, linewidth=1.5, label=r"$x(\beta)$")
    ax1.set_ylabel(r"$x$ (graus)", fontsize=12)
    ax1.set_title(r"Solução $x$ em função de $\beta$", fontsize=12)
    ax1.legend()
    ax1.grid(True, linestyle="--", alpha=0.6)

    ax2.bar(betas_deg, iteracoes, width=(betas_deg[1] - betas_deg[0]) * 0.8,
            color="#10B981", edgecolor="#065F46", alpha=0.85, label="Iterações")
    ax2.set_xlabel(r"$\beta$ (graus)", fontsize=12)
    ax2.set_ylabel("Nº de iterações", fontsize=12)
    ax2.set_title("Iterações por valor de β", fontsize=12)
    ax2.legend()
    ax2.grid(True, linestyle="--", alpha=0.6, axis="y")

    plt.tight_layout()
    plt.savefig("freudenstein_resultados.png", dpi=150, bbox_inches="tight")
    plt.show()
    print("\n✔ Gráfico salvo em 'freudenstein_resultados.png'")


def main():
    print(f"l1={L1}, l2={L2}, l3={L3}, l4={L4} cm")
    print(f"β ∈ [0°, {np.degrees(BETA_MAX):.1f}°] | {NUM_BETA} valores | tol={TOL_REL:.0e} | max_iter={MAX_ITER}")

    betas = np.linspace(BETA_MIN, BETA_MAX, NUM_BETA)
    solucoes, iteracoes, convergencias = resolver(betas)

    # Tabela de resultados
    print(f"\n{'β (°)':>10} {'x (°)':>10} {'x (rad)':>12} {'Iter':>6} {'Conv':>5}")
    print("-" * 50)
    for b, s, it, c in zip(betas, solucoes, iteracoes, convergencias):
        print(f"{np.degrees(b):10.3f} {np.degrees(s):10.3f} {s:12.6f} {it:6d} {'✔' if c else '✘':>5}")

    print(f"\nConvergiram: {np.sum(convergencias)}/{len(betas)}")

    plotar(betas, solucoes, iteracoes)


if __name__ == "__main__":
    main()
