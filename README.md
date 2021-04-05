# MAP3122-EP02

## Repositório criado para resolver o Exercício Programa 2 da disciplina MAP3122 - Métodos Numéricos e Aplicações

Os programas dos exercícios computacionais foram escritos em python 3, na versão 3.8.5, e usando o pacote numpy. A entrada e saída foram feitas de forma a ajudar o usuário a executar o programa e facilitar a análise dos resultados. E a biblioteca matplotlib foi usada para as plotagens.

Para a realização do exercício foram criados 7 arquivos `.py` os arquivos `exercicio1.py`, `exercicio2.py` e `exercicio3.py`, sendo cada um respectivo a um dos exercícios. Além disso, foram criados os arquivos `rk4.py`, `euler_forward.py`, `euler_backward.py` e `plotter.py` contendo a implementação do seu respectivo método de resolução de EDO e a última contém métodos de plotagem de gráficos comuns aos exercícios. 

Para rodar o exercício programa, deve-se ter python3 instalado no computador onde se irá executar. Além disso, os scripts possuem dependências, as quais estão listadas no arquivo `requirements.txt` e podem ser instaladas ao se rodar no terminal:

```bash
pip3 install -r requirements.txt

```

Por fim, para se rodar o exercício programa, estando no mesmo diretório do arquivo (pasta `src/`), deve-se rodar no terminal:

```bash
python3 [NOME DO ARQUIVO DO PROGRAMA]

```

Você deve substituir [NOME DO ARQUIVO DO PROGRAMA] por `exercicio1.py` ou `exercicio2.py` ou `exercicio3.py`, de acordo com o exercício que se deseja rodar.

Então para rodar o programa exercicio1.py, por exemplo, basta fazer:
```bash
python3 exercicio1.py

```