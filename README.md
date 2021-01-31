# genomics
Projeto para aprendizado da disciplina de Ciência de Dados

## Obtendo o projeto
### Usando Git
* Instale o git em sua máquina
* Abra o terminal ou prompt na pasta onde deseja baixar o projeto
* Execute o comando:
```shell
git clone https://github.com/SousaPedro11/genomics.git
```

### Usando a interface Web
* Clique no botão verde escrito CODE
* Clique na opção Download ZIP
* Após o Download extraia o conteúdo do arquivo .zip

## Preparação do ambiente
* Pelo terminal ou Prompt vá até a pasta raiz do projeto
* Instale o Python 3
* Instale o Anaconda ou miniconda
* Crie o ambiente virtual a partir de [environment.yml](environment.yml):
```shell
conda env create --name env --file=environments.yml
```

## Executando o Projeto
* Ative o ambiente virtual:
```shell
conda activate ./env
```
* Execute o projeto a partir do arquivo [main.py](main.py):
```shell
python main.py
```
