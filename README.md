# Executando o Solver da Equação de Freudenstein

Instruções diretas para configurar o ambiente, isolar as dependências e executar o código.

## 1. Definir Versão do Python

Defina a versão do Python local com `pyenv`:

```bash
# Instalar a versão caso não tenha
pyenv install 3.10.12

# Definir a versão local do projeto
pyenv local 3.10.12
```

## 2. Isolar e Instalar Dependências

Crie um ambiente virtual (`venv`) local para garantir que as dependências fiquem isoladas do sistema:

```bash
# Criar o ambiente virtual na pasta do projeto
python -m venv .venv

# Ativar o ambiente virtual
source .venv/bin/activate

# Instalar as dependências isoladamente (incluindo pygame para a animação)
pip install numpy matplotlib pygame
```

## 3. Executar o Código

Com o ambiente virtual ativo, você pode executar tanto o solver numérico quanto a animação interativa.

### Executar o Solver e Gerar Dados/Tabela
Para rodar o script principal que gera o gráfico estático e a tabela de convergência em LaTeX:
```bash
python freudenstein_newton_raphson.py
```

> [!IMPORTANT]
> ### 🎮 Executar a Animação Interativa (Pygame)
> Para visualizar o movimento do mecanismo de 4 barras em tempo real e interagir com o ângulo de entrada ($\beta$) e a aproximação inicial ($x_0$), execute o script de animação:
> 
> ```bash
> python anim.py
> ```
> 
> **Funcionalidades da Animação:**
> - **Controle de Sliders:** Ajuste dinamicamente o ângulo de entrada $\beta$ e o chute inicial $x_0$ usando a barra de controle na parte inferior.
> - **Diagnóstico em Tempo Real:** O painel informa se o solver convergiu, o número de iterações necessárias e o ângulo de saída $\alpha$.
> - **Feedback Visual de Convergência:** O mecanismo é desenhado na tela somente se o solver convergir. Caso contrário, um aviso visual vermelho de falha é exibido.

*Para desativar o ambiente virtual quando terminar, execute:* `deactivate`

## 4. Inserindo a Tabela no Overleaf

A execução gera o arquivo `freudenstein_tabela.tex` com a tabela de resultados formatada em LaTeX.

**Pacotes necessários** — adicione ao preâmbulo do seu documento:

```latex
\usepackage{longtable}  % tabelas multi-página
\usepackage{booktabs}   % réguas horizontais profissionais (\toprule, \midrule, \bottomrule)
```

**Inserindo a tabela** — no corpo do documento, onde desejar:

```latex
\input{freudenstein_tabela.tex}
```

Referência cruzada: `Tabela~\ref{tab:freudenstein}`