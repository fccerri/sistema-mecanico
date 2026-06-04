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

# Instalar as dependências isoladamente
pip install numpy matplotlib
```

## 3. Executar o Código

Com o ambiente virtual ativo, execute o script principal:

```bash
python freudenstein_newton_raphson.py
```

*Para desativar o ambiente virtual quando terminar, execute:* `deactivate`